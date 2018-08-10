#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "learn.h"
#include "matrix.h"

int extractMatRow(Matrix *extract, Matrix *mat, int row){
    if (row > (mat->rows)-1){
        return -1;
    }

    int rows = 1;
    int cols = mat->cols;

    emptyMat(extract, rows, cols);

    int c = 0;
    while (c < cols){
        extract->ops[0][c] = mat->ops[row][c];
        c++;
    }

    return 0;
}

int concatOnes(Matrix *concat, Matrix *mat){
    Matrix *ones = (Matrix*) malloc(sizeof(Matrix));
    emptyMat(ones, mat->rows, 1);
    int r = 0;
    while (r < ones->rows){
        ones->ops[r][0] = 1;
        r++;
    }

    concatMatLR(concat, ones, mat);
    
    freeMatStruct(ones);
    free(ones);

    return 0;
}

void extractTestMeta(char *srcStr, int *points, char **waste){
    char *lnTkn = strtok_r(srcStr, "\n", &srcStr);
    
    *points = atoi(lnTkn);

    int wasteLen = strlen(srcStr);
    waste[0] = (char*)malloc( (wasteLen + 1) * sizeof(char) );
    waste[0][0] = '\0';
    strncat(waste[0], srcStr, wasteLen);
}

void extractTrainMeta(char *srcStr, int *attr, int *exmpl, char **waste){
    char *lnTkn = strtok_r(srcStr, "\n", &srcStr);
    
    *attr = atoi(lnTkn);
    lnTkn = strtok_r(NULL, "\n", &srcStr);
    *exmpl = atoi(lnTkn);

    int wasteLen = strlen(srcStr);
    waste[0] = (char*)malloc( (wasteLen + 1) * sizeof(char) );
    waste[0][0] = '\0';
    strncat(waste[0], srcStr, wasteLen);
}

char *readFile(char *filename){
    FILE *file = NULL;
    file = fopen(filename, "r");
  
    if (file == NULL){
        //printf("error\n");
        return NULL; // failed to open file
    }
  
    int fileContSize = 1;
    char *fileContent = (char*) malloc(fileContSize * sizeof(char));
    if (fileContent == NULL){
        free(fileContent);
        fclose(file);
        return NULL; // failed to allocate memory
    }
    fileContent[0] = '\0';
    char buffer[1024];
    buffer[0] = '\0';
  
    while ( fgets(buffer, 1024, file) != NULL ){
        int contSizeDiff = fileContSize - (fileContSize + (strlen(buffer) + 1));
        if (contSizeDiff < 0){
            contSizeDiff *= -1;
            fileContent = realloc(fileContent, fileContSize + contSizeDiff);
            fileContSize += contSizeDiff;
            if (fileContent == NULL){
                free(fileContent);
                fclose(file);
                return NULL; // failed to allocate memory
            }
        }
        strncat(fileContent, buffer, 1023);
    }
  
    fclose(file);
  
    //printf("Read In: %s\n", fileContent);
  
    return fileContent;
}
  

int main(int argc, const char **argv){
    
    // --- READ FILENAME --- //
    char *trainFname = (char*)malloc( (strlen(argv[1]) + 1) * sizeof(char) );
    if (trainFname == NULL){
        return 0; // failed to read cmd args
    }
    trainFname[0] = '\0';
    snprintf(trainFname, strlen(argv[1]) + 1, "%s", argv[1]);

    char *testFname = (char*)malloc( (strlen(argv[2]) + 1) * sizeof(char) );
    if (testFname == NULL){
        free(trainFname);
        return 0; // failed to read cmd args
    }
    testFname[0] = '\0';
    snprintf(testFname, strlen(argv[2]) + 1, "%s", argv[2]);

    
    // --- OPEN AND READ FILE --- //
    char *trainFC = readFile(trainFname);
    if (trainFC == NULL){
        free(trainFname);
        free(testFname);
        return 0; // failed to read file
    }
    free(trainFname);
    
    char *testFC = readFile(testFname);
    if (testFC == NULL){
        free(testFname);
        free(trainFC);
        return 0; // failed to read file
    }
    free(testFname);

    // --- GRAB METADATA (TRAIN) --- //
    int K = 0; // attributes
    int N = 0; // number of training examples
    char *matrixStr_train;
    extractTrainMeta(trainFC, &K, &N, &matrixStr_train);

    // --- BUILD MATRIX FROM FILE --- //
    Matrix *f = (Matrix*) malloc(sizeof(Matrix));
    buildMat(matrixStr_train, f);

    free(matrixStr_train);
    free(trainFC);

    // --- GRAB METADATA (TEST) --- //
    int M = 0; // number of data points
    char *matrixStr_test;
    extractTestMeta(testFC, &M, &matrixStr_test);

    // --- BUILD MATRIX FROM FILE --- //
    Matrix *f_test = (Matrix*) malloc(sizeof(Matrix));
    buildMat(matrixStr_test, f_test);

    free(matrixStr_test);
    free(testFC);
    
    // --- CONCATENATE ONES COLUMN TO FIRST COLUMN --- //

    Matrix *con_train = (Matrix*) malloc(sizeof(Matrix));
    concatOnes(con_train, f);
    
    freeMatStruct(f);
    free(f);

    Matrix *con_test = (Matrix*) malloc(sizeof(Matrix));
    concatOnes(con_test, f_test);
    
    freeMatStruct(f_test);
    free(f_test);

    // --- EXTRACT Y COLUMN (from TRAIN) --- //
    Matrix *Y = (Matrix*) malloc(sizeof(Matrix));
    emptyMat(Y, con_train->rows, 1);
    int r = 0;
    while (r < Y->rows){
        Y->ops[r][0] = con_train->ops[r][(con_train->cols)-1];
        r++;
    }

    Matrix *X = (Matrix*) malloc(sizeof(Matrix));
    delMatCol(X, con_train, (con_train->cols)-1);

    freeMatStruct(con_train);
    free(con_train);

    // --- CALCULATE W --- //
    // W = C * Y
    // C = B * trans(X)
    // B = inv(A)
    // A = trans(X) * X

    Matrix *W = (Matrix*) malloc(sizeof(Matrix));
    Matrix *C = (Matrix*) malloc(sizeof(Matrix));
    Matrix *B = (Matrix*) malloc(sizeof(Matrix));
    Matrix *A = (Matrix*) malloc(sizeof(Matrix));
    Matrix *t_X = (Matrix*) malloc(sizeof(Matrix));

    //puts("X MATRIX: ");
    //printMat(X);
    
    transposeMat(t_X, X);
    //puts("TRANSPOSE OF X: ");
    //printMat(t_X);

    multiplyMat(A, t_X, X);
    //puts("A: ");
    //printMat(A);

    inverseMat(B, A);
    //puts("B: ");
    //printMat(B);

    multiplyMat(C, B, t_X);
    //puts("C: ");
    //printMat(C);

    //puts("Y MATRIX: ");
    //printMat(Y);

    multiplyMat(W, C, Y);
    //puts("W: ");
    //printMat(W);
    
    //puts("DONE!");

    freeMatStruct(C);
    freeMatStruct(B);
    freeMatStruct(A);
    freeMatStruct(t_X);
    freeMatStruct(X);
    freeMatStruct(Y);
    free(C);
    free(B);
    free(A);
    free(t_X);
    free(X);
    free(Y);

    // --- CALCULATE PREDICTION PER ROW --- //
    Matrix *Y_predict = (Matrix*) malloc(sizeof(Matrix));
    emptyMat(Y_predict, con_test->rows, 1);
    r = 0;
    while (r < con_test->rows){
        Matrix *extr = (Matrix*) malloc(sizeof(Matrix));
        extractMatRow(extr, con_test, r);
        Matrix *pred = (Matrix*) malloc(sizeof(Matrix));
        multiplyMat(pred, extr, W);
        Y_predict->ops[r][0] = pred->ops[0][0];
        freeMatStruct(extr);
        free(extr);
        freeMatStruct(pred);
        free(pred);
        r++;
    }

    //puts("TEST X: ");
    //printMat(con_test);

    //puts("Y PREDICT: ");
    printPredictMat(Y_predict);

    freeMatStruct(W);
    free(W);

    freeMatStruct(con_test);
    free(con_test);

    freeMatStruct(Y_predict);
    free(Y_predict);

    return 0;
}

