// Para compilar con mpi: mpic++ main.cpp
// Correr con: time mpirun -np 2 a.out

#include <iostream>
#include <fstream>
using namespace std;
#include <mpi.h>
#include <cmath>

// Funcion para leer dimensiones de matriz.
void dims(const string &fname, int &nrows, int &ncols){
    // ifstream file(fname);
    ifstream file(fname.c_str());
    if (file.is_open()) {
        file >> nrows;
        file >> ncols;
        file.close();
    }
}

// Funcion que calcula la cantidad de elementos que cada segmento (de Matriz o Vector) va a tener.
void len_segmento(int &len_segmento_vector, int &len_segmento_matriz, int nrows, int b_len, int world_size, int world_rank){
    // Cantidad de elementos de Vector:
    int cantidad_vector = b_len / world_size;
    int indice_rango_vector = world_rank * cantidad_vector;
    if(world_rank == world_size-1){
        cantidad_vector += b_len % world_size;
    }
    len_segmento_vector = cantidad_vector;

    // Cantidad de elementos de Matriz:
    int cantidad_matriz = nrows / world_size;
    int indice_rango_matriz = world_rank * cantidad_matriz;
    if(world_rank == world_size-1){
        cantidad_matriz += nrows % world_size;
    }
    len_segmento_matriz = cantidad_matriz;
}

// Funcion para leer segmento de .txt.
double* leer_segmento(const string &fname, int inicio, int fin) {
    // ifstream file(fname);
    ifstream file(fname.c_str()); 
    // double* array = nullptr;
    double* array = NULL;

    if (file.is_open()) {
        int nrows, ncols;
        file >> nrows;
        file >> ncols;

        for(int i = 0; i < inicio; i++) {
            double tmp;
            file >> tmp;
        }
        
        int intervalo = (fin - inicio);
        // array = (double*)malloc(intervalo * sizeof(double));
        array = new double[intervalo];

        // if (array != nullptr) {
        if (array != NULL) {
            for (int i = 0; i < intervalo; i++) {
                file >> array[i];
            }
        }
        file.close();
    }

    return array;
}


void indices(int &inicio, int &fin, int world_size, int world_rank, int ncols, int nrows, int len_segmento_matriz){
    if(world_rank == 0){
        inicio = 0;
        fin = ncols*len_segmento_matriz;
    }
    else if(world_rank == world_size-1){
        inicio = ncols*(nrows/world_size)*(world_size-1);
        fin = inicio+ncols*len_segmento_matriz;
    }
    else{
        inicio = ncols*len_segmento_matriz + (ncols*len_segmento_matriz*(world_rank-1));
        fin = inicio+ncols*len_segmento_matriz;
    }
}

// Funcion que calcula valor de numerador de vector b_(k+1)
double* b_k_1(double* array_valores, double* final_array, int len_segmento_matriz, int b_len, int world_rank){
    double* b_k_mas_1 = new double[len_segmento_matriz];
    int indice = 0;
    for(int i = 0; i < len_segmento_matriz; i++){
        double suma_fila = 0.0;
        for(int j = 0; j < b_len; j++){
            suma_fila += array_valores[indice] * final_array[j];
            indice += 1;
        }
        b_k_mas_1[i] = suma_fila;
    }
    return b_k_mas_1;
}

// Funcion que calcula vector b_(k+1) multiplicado por el denominador del enunciado.
double* b_k_1_dist(double* b_k_mas_1, int len_segmento_matriz){
    double suma_de_cuadrados = 0.0;
    for(int i = 0; i < len_segmento_matriz; i++){
        suma_de_cuadrados += b_k_mas_1[i] * b_k_mas_1[i];
    }
    double norma_2 = 1.0 / std::sqrt(suma_de_cuadrados);

    for(int i = 0; i < len_segmento_matriz; i++){
        b_k_mas_1[i] *= norma_2;
    }
    return b_k_mas_1;
}

// Funcion que calcula valor propio mas grande final.
double u_k(double* b_k_mas_1_temp, double* matriz, double* v_transpuesto, int b_len, int nrows, int ncols){
    // denominador
    double denominador = 0;
    for(int i = 0; i < b_len; i++){
        denominador += v_transpuesto[i] * b_k_mas_1_temp[i];
    }
    denominador = 1/denominador;

    // numerador
    // double* multi_t_matriz = new double[nrows];
    // double* multi_t_matriz = (double*)malloc(nrows * sizeof(double));
    double* multi_t_matriz = new double[nrows];
    int indice = 0;
    for(int i = 0; i < nrows; i++){
        double temp = 0.0;
        for(int j = 0; j < ncols; j++){
            temp += v_transpuesto[j] * matriz[indice];
            indice += 1;
        }
        multi_t_matriz[i] = temp;
    }

    double valor_final = 0.0;
    for(int i=0; i < nrows; i++){
        valor_final += multi_t_matriz[i] * b_k_mas_1_temp[i];
    }

    valor_final *= denominador;
    return valor_final;
}

// Funcion que consigue transpuesto de vector.
double* transpuesto(double* array_valores, int b_len){
    double* v_transpuesto = new double[b_len];
    // double* v_transpuesto = (double*)malloc(b_len * sizeof(double));
    int c = 0;
    for(int i = b_len-1; i >= 0; i--){
        v_transpuesto[c] = array_valores[i];
        c += 1;
    }
    return v_transpuesto;
}

int main()
{
    // Conseguimos dimensiones matriz y longitud vector b_0.
    int nrows, ncols, b_len;
    dims("matrix_1000.txt", nrows, ncols);
    b_len = ncols;

    // Inicializamos MPI y conseguimos cantidad procesos y rango de cada uno.
    MPI_Init(NULL,NULL);
    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Imprimimos nombre de procesador en cada proceso.
    int err;
    char procname[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    err = MPI_Get_processor_name(procname, &name_len);
    cout << "Process number: " << world_rank << " Processor name: " << procname << endl;

    // Obtenemos cantidad de filas que cada proceso debe tener, y cantidad elementos que cada proceso debe generar del vector
    int len_segmento_vector, len_segmento_matriz;
    len_segmento(len_segmento_vector, len_segmento_matriz, nrows, b_len, world_size, world_rank);

    // Cada proceso genera su parte del vector
    double local_vector[len_segmento_vector]; // Buffer inicial
    for(int x = 0; x < len_segmento_vector; x++){
        local_vector[x] = 1.0;
    }

    // All Gather no me funciono, tengo entendido que es porque no siempre la longitud de los vectores a 'traspasar' son iguales.
    // Debido a esto, hay que ocupar Allgatherv, aqui la documentacion: https://www.open-mpi.org/doc/v4.1/man3/MPI_Allgatherv.3.php

    double* final_array = new double[b_len];
    int recvcounts[world_size];
    int displs[world_size];

    for (int i = 0; i < world_size; i++) {
        if(i == world_size-1){
            recvcounts[i] = (b_len/world_size) + (b_len%world_size);
            displs[i] = (b_len/world_size)*i;
        }
        else{
            recvcounts[i] = (b_len/world_size);
            displs[i] = i * recvcounts[i];
        }
    }
    MPI_Allgatherv(local_vector, len_segmento_vector, MPI_DOUBLE, final_array, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD);

    // Indices para seleccionar filas.
    int inicio, fin;
    indices(inicio, fin, world_size, world_rank, ncols, nrows, len_segmento_matriz);

    for(int i = 0; i < b_len; i++){
        printf("P%d, indice: %d, valor: %f\n", world_rank, i, final_array[i]);
    }

    // // AQUI
    // // Array 1D que contiene valores de fila/s matriz.
    // double* array_valores = leer_segmento("matrix_1000.txt", inicio, fin);

    // // Funcion que calcula valor de b_(k+1), luego iteramos por una cantidad 'n'.
    // double* init_buff = new double[nrows];
    // double* b_k_mas_1_num = b_k_1(array_valores, final_array, len_segmento_matriz, b_len, world_rank);
    // MPI_Allgatherv(b_k_mas_1_num, len_segmento_matriz, MPI_DOUBLE, init_buff, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD);
    // double* b_k_mas_1 =  b_k_1_dist(init_buff, nrows);

    // // Iteramos y actualizamos el vector para un numero 'n' fijo de iteraciones.
    // int n = 200;
    // double* b_k_mas_1_temp = new double[b_len];
    // for(int i = 0; i < n; i++){
    //     b_k_mas_1 =  b_k_1(array_valores, b_k_mas_1, len_segmento_matriz, b_len, world_rank);
    //     MPI_Allgatherv(b_k_mas_1, len_segmento_matriz, MPI_DOUBLE, b_k_mas_1_temp, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD);
    //     b_k_mas_1 =  b_k_1_dist(b_k_mas_1_temp, nrows);
    // }

    // // Finalmente con el vector propio ya "convergido", Solamente proceso 0 lee toda la matriz y 
    // // la multiplica por valor final de vector b_(k+1) y su transpuesto.
    // if(world_rank == 0){
    //     double* matriz = leer_segmento("matrix_1000.txt", 0, (nrows*ncols));
    //     double* v_transpuesto = transpuesto(b_k_mas_1_temp, b_len);
    //     double v_propio = u_k(b_k_mas_1_temp, matriz, v_transpuesto, b_len, nrows, ncols);
    //     double valor_final = 10 - v_propio;
    //     printf("VALOR FINAL: %f\n", valor_final);
    // }

    MPI_Finalize();

    return 0;
}