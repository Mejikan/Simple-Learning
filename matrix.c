#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "matrix.h"

int delMatCol(Matrix *res, Matrix *src, int col){
    if (col >= (src->cols)){
        return -1;
    }
    int rows = src->rows;
    int cols = (src->cols)-1;
    
    double **matrix = (double**)malloc(rows * sizeof(double*)); // reserve space for row pointers
    if (matrix == NULL){
        return -1;
    }
    
    int r = 0;
    while (r < rows){
        matrix[r] = (double*)malloc(cols * sizeof(double)); //reserve space for cols
        if (matrix[r] == NULL){ // failed to allocate memory
            freeMat(matrix, r);
            return -1;
        }
        r++;
    }

    r = 0;
    while (r < rows){
        int c = 0;
        int ci = 0;
        while (c < (src->cols)){
            if (c != col){
                matrix[r][ci] = src->ops[r][c];
                ci++;
            }
            c++;
        }
        r++;
    }
    
    
    res->rows = rows;
    res->cols = cols;
    res->ops = matrix;

    return 0;
}

int multiplyMat(Matrix *prod, Matrix *matL, Matrix *matR){
    if (matL->cols != matR->rows){
        return -1;
    }

    int rows = matL->rows;
    int cols = matR->cols;

    double **matrix = (double**)malloc(rows * sizeof(double*)); // reserve space for row pointers
    if (matrix == NULL){
        return -1;
    }

    int r = 0;
    while (r < rows){
        matrix[r] = (double*)malloc(cols * sizeof(double)); //reserve space for cols
        if (matrix[r] == NULL){ // failed to allocate memory
            freeMat(matrix, r);
            return -1;
        }
        r++;
    }

    r = 0;
    while(r < rows){
        int c = 0;
        while(c < cols){

            double sum = 0;
            int ci = 0; // column iterator
            while (ci < matL->cols){
                double prod = (matL->ops[r][ci]) * (matR->ops[ci][c]);
                sum += prod;
                ci++;
            }

            matrix[r][c] = sum;
            c++;
        }
        r++;
    }

    prod->rows = rows;
    prod->cols = cols;
    prod->ops = matrix;

    return 0;
}

int transposeMat(Matrix *trans, Matrix *mat){
    int rows = mat->cols;
    int cols = mat->rows;
    double **matrix = (double**)malloc(rows * sizeof(double*)); // reserve space for row pointers
    if (matrix == NULL){
        return -1;
    }

    int r = 0;
    while (r < rows){
        matrix[r] = (double*)malloc(cols * sizeof(double)); //reserve space for cols
        if (matrix[r] == NULL){ // failed to allocate memory
            freeMat(matrix, r);
            return -1;
        }

        // insert values
        int c = 0;
        while (c < cols){
            matrix[r][c] = mat->ops[c][r];
            c++;
        }
        r++;
    }

    trans->rows = rows;
    trans->cols = cols;
    trans->ops = matrix;

    return 0;
}

int inverseMat(Matrix *inv, Matrix *mat){
    // use Gauss-Jordan
    Matrix *ident = (Matrix*) malloc(sizeof(Matrix));
    if (ident == NULL){
        return -1;
    }
    Matrix *con = (Matrix*) malloc(sizeof(Matrix));
    if (con == NULL){
        free(ident);
        return -1;
    }

    if (identMat(ident, mat->rows) != 0){
        free(ident);
        return -1;    
    }

    if (concatMatLR(con, mat, ident) != 0){
        free(ident);
        free(con);
        return -1;    
    }

    // elimination process
    int matCols = mat->cols;
    int r = 0;
    while (r < con->rows){
        int c = 0;
        while (c < matCols){
            double *ei = &con->ops[r][(matCols)+c]; //pointer to element in identity matrix
            double eiv = (*ei); // value of ei
            double *em = &con->ops[r][c]; // pointer to element in original matrix
            double emv = (*em); // value of em
            if ( r == c ){
                // divide by emv
                int ci = c; // column iterator index for entire row
                while (ci < con->cols){
                    con->ops[r][ci] = (con->ops[r][ci])/emv;
                    ci++;
                }
                break; // stop column traversing here
            } else {
                // subtract by (emv * pivot row)
                // the column currently in has a pivot that exists at the row equal to column number
                double *pivotRow = con->ops[c];

                int ci = c; // column iterator index for entire row
                while (ci < con->cols){
                    con->ops[r][ci] = (con->ops[r][ci]) - (emv * pivotRow[ci]);
                    ci++;
                }
            }
            c++;
        }
        r++;
    }

    r = r-1;
    while (r >= 0){
        int c = matCols-1;
        while (c >= 0){
            double *ei = &con->ops[r][(matCols)+c]; //pointer to element in identity matrix
            double eiv = (*ei); // value of ei
            double *em = &con->ops[r][c]; // pointer to element in original matrix
            double emv = (*em); // value of em
            if (r == c){
                break; // stop column traversing here
            } else {
                // subtract by (emv * pivot row)
                // the column currently in has a pivot that exists at the row equal to column number
                double *pivotRow = con->ops[c];

                int ci = c; // column iterator index for entire row
                while (ci < con->cols){
                    con->ops[r][ci] = (con->ops[r][ci]) - (emv * pivotRow[ci]);
                    ci++;
                }
            }
            c--;
        }
        r--;
    }
    
    r = 0;
    while (r < ident->rows){
        int c = 0;
        while (c < ident->cols){
            ident->ops[r][c] = con->ops[r][matCols+c];
            c++;
        }
        r++;
    }
    
    inv->rows = ident->rows;
    inv->cols = ident->cols;
    inv->ops = ident->ops;
    
    freeMatStruct(con);
    free(con);
    free(ident);

    return 0;
}

int concatMatLR(Matrix *concat, Matrix *matL, Matrix *matR){
    if (matL->rows != matR->rows){
        return -1;
    }

    int rows = matL->rows;
    int cols = (matL->cols) + (matR->cols);
    double **matrix = (double**)malloc(rows * sizeof(double*)); // reserve space for row pointers
    if (matrix == NULL){
        return -1;
    }

    int r = 0;
    while (r < rows){
        matrix[r] = (double*)malloc(cols * sizeof(double)); //reserve space for cols
        if (matrix == NULL){ // failed to allocate memory
            freeMat(matrix, rows);
            return -1;
        }

        // insert values
        int c = 0;
        int c1 = 0;
        while (c1 < matL->cols){
            matrix[r][c] = matL->ops[r][c1];
            c++;
            c1++;
        }
        c1 = 0;
        while (c1 < matR->cols){
            matrix[r][c] = matR->ops[r][c1];
            c++;
            c1++;
        }
        r++;
    }

    concat->rows = rows;
    concat->cols = cols;
    concat->ops = matrix;

    return 0;
}

void printMat(Matrix *mat){
    int r = 0;
    while (r < mat->rows){
        int c = 0;
        while (c < mat->cols){
            printf("%lf ", mat->ops[r][c]);
            c++;
        }
        printf("\n");
        r++;
    }
}

void printPredictMat(Matrix *mat){
    int r = 0;
    while (r < mat->rows){
        int c = 0;
        while (c < mat->cols){
            printf("%0.0f\n", mat->ops[r][c]);
            c++;
        }
        r++;
    }
}

/**
 * Builds the identity of a square matrix
 * and inserts it into the matrix given.
 **/
int identMat(Matrix *mat, int len){
    int rows = len;
    int cols = len;
    double **matrix = (double**)malloc(rows * sizeof(double*)); // reserve space for row pointers
    if (matrix == NULL){ // failed to allocate memory
        return -1;
    }

    int r = 0;
    while (r < rows){
        matrix[r] = (double*)malloc(cols * sizeof(double)); //reserve space for cols
        if (matrix == NULL){ // failed to allocate memory
            freeMat(matrix, rows);
            return -1;
        }

        // insert values
        int c = 0;
        while (c < cols){
            if (r == c){
                matrix[r][c] = 1;
            } else {
                matrix[r][c] = 0;
            }
            c++;
        }
        r++;
    }

    mat->rows = rows;
    mat->cols = cols;
    mat->ops = matrix;

    return 0;
}

void freeMatStruct(Matrix *mat){
    freeMat(mat->ops, mat->rows);
}

void freeMat(double **mat, int rows){
    int r = 0;
    while( r < (rows) ){
        free(mat[r]);
        r++;
    }
    free(mat);
}

int buildMat(char *srcStr, Matrix *mat){
    char *rowTkn = strtok_r(srcStr, "\n", &srcStr);
    
    double **matrix = (double**)malloc(1 * sizeof(double*)); // reserve space for 1 int ptr: 1 row, 1 col
    if (matrix == NULL){ // failed to allocate memory
        return -1;
    }
    
    int rows = 0; // number of total rows
    int cols = 0; // number of total columns
    
    int rowIdx = 0; // starts at index 1 for each row - counter
    
    while (rowTkn != NULL){
        matrix = (double**)realloc(matrix, (rowIdx+1) * sizeof(double*));
        if (matrix == NULL){ // failed to allocate memory
            freeMat(matrix, rowIdx);
            return -1;
        }

        if (rowIdx == 0){ // first row
            matrix[rowIdx] = (double*)malloc(1 * sizeof(double)); // allocate space for 1 * sizeof double, for 1 column
            char *colTkn = strtok_r(rowTkn, ",", &rowTkn);
            int colIdx = 0;
            while (colTkn != NULL){
                matrix[rowIdx] = (double*)realloc(matrix[rowIdx], (colIdx+1) * sizeof(double));
                if (matrix[rowIdx] == NULL){ // failed to allocate memory
                    freeMat(matrix, rowIdx);
                    return -1;
                }
                matrix[rowIdx][colIdx] = atof(colTkn);
                colIdx++;
                colTkn = strtok_r(NULL, ",", &rowTkn);
            }
            cols = colIdx;
        } else {
            matrix[rowIdx] = (double*)malloc(cols * sizeof(double)); // allocate space for cols * sizeof double
            if (matrix[rowIdx] == NULL){ // failed to allocate memory
                freeMat(matrix, rowIdx);
                return -1;
            }
            char *colTkn = strtok_r(rowTkn, ",", &rowTkn);
            int colIdx = 0;
            while (colTkn != NULL && colIdx < cols){
                matrix[rowIdx][colIdx] = atof(colTkn);
                colIdx++;
                colTkn = strtok_r(NULL, ",", &rowTkn);
            }
        }
        
        rowIdx++;
        rowTkn = strtok_r(NULL, "\n", &srcStr);
    }

    rows = rowIdx;

    mat->rows = rows;
    mat->cols = cols;
    mat->ops = matrix;

    return 0;
}

int emptyMat(Matrix *mat, int rows, int cols){
    double **matrix = (double**)malloc(rows * sizeof(double*)); // reserve space for row pointers
    if (matrix == NULL){
        return -1;
    }

    int r = 0;
    while(r < rows){
        matrix[r] = (double*)malloc(cols * sizeof(double)); // allocate space for cols * sizeof double
        if (matrix[r] == NULL){
            freeMat(matrix, r);
            return -1;
        }
        r++;
    }

    mat->rows = rows;
    mat->cols = cols;
    mat->ops = matrix;

    return 0;
}
