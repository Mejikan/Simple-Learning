typedef
    struct {
        int rows;
        int cols;
        double **ops;
    }
Matrix;

int delMatCol(Matrix*, Matrix*, int);

int multiplyMat(Matrix*, Matrix*, Matrix*);

int transposeMat(Matrix*, Matrix*);

int inverseMat(Matrix*, Matrix*);

int concatMatLR(Matrix*, Matrix*, Matrix*);

void printMat(Matrix*);

void printPredictMat(Matrix*);

int identMat(Matrix*, int);

void freeMatStruct(Matrix*);

void freeMat(double**, int);

int buildMat(char*, Matrix*);

int emptyMat(Matrix*, int, int);