/* linalg_lapack_nif.c */

#include <stdlib.h>

#include "cblas.h"
#include "erl_nif.h"

extern void dgesdd_(const char *jobz,
				const int *m, const int *n,
				double *a, const int *lda, double *s,
				double *u, const int *ldu,
				double *vt, const int *ldvt,
				double *work, const int *lwork, int *iwork, int *info);
 

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

static ERL_NIF_TERM transpose_nif(ErlNifEnv* env, int argc, const ERL_NIF_TERM argv[])
{

   unsigned int ncolA=-1,nrowA=-1;
   double val;

    ERL_NIF_TERM list, head, tail, sublist, subhead,subtail;
    list=argv[0];

    enif_get_list_cell(env, list, &head, &tail);
    if (!enif_get_list_length(env, list, &nrowA))
        return enif_make_badarg(env);
    if(!enif_get_list_length(env, head, &ncolA))
        return enif_make_badarg(env);

    if (nrowA==0)
        return enif_make_badarg(env);

    double *a = enif_alloc(ncolA*nrowA * sizeof *a);
    ERL_NIF_TERM *A = enif_alloc(ncolA*nrowA * sizeof *A);

    double *pcoeffa;
    pcoeffa=&a[0];

    ERL_NIF_TERM *pcoeffA;
    pcoeffA=&A[0];

    while (enif_get_list_cell(env, list, &head, &tail)) {

        sublist=head;
        while (enif_get_list_cell(env, sublist, &subhead, &subtail)) {
            if (!enif_get_double(env, subhead, &val)) 
				return enif_make_badarg(env);

            *pcoeffa = val;
            pcoeffa++;

            *pcoeffA =enif_make_double(env,val);
            pcoeffA++;

            sublist = subtail;
        }
        list = tail;
    }
    *pcoeffa = '\0';
    *pcoeffA = '\0';

	return enif_make_list_from_array(env, A, ncolA*nrowA);

/*
    list=argv[1];

    int i,j;

    int nrowC=nrowA,ncolC=ncolB;
    double *c = enif_alloc(ncolB*nrowA * sizeof *c);
    for (i=0; i<nrowC; i++) {
    	for (j=0; j<ncolC; j++) {
		 	c[(i*ncolA)+j]=0;
		}
 	}

    ERL_NIF_TERM *C = enif_alloc(ncolB*nrowA * sizeof *C);
	ERL_NIF_TERM c_nif = enif_make_list(env, 0);

    for (i=nrowC-1; i>=0; i--) {
    	for (j=0; j<ncolC; j++) {
		 	C[j]=enif_make_double(env,c[(i*ncolC)+j]);
		}
		c_nif=enif_make_list_cell(env,enif_make_list_from_array(env, C, ncolC),c_nif);
 	}
*/

    // and then clean them up too:
//    enif_free(a);
//    enif_free(A);

//	return c_nif;
}

static ERL_NIF_TERM matmul_nif(ErlNifEnv* env, int argc, const ERL_NIF_TERM argv[])
{

   unsigned int ncolA=-1,nrowA=-1,ncolB=-1,nrowB=-1;
   double val;

    ERL_NIF_TERM list, head, tail, sublist, subhead,subtail;
    list=argv[0];

    enif_get_list_cell(env, list, &head, &tail);
    if (!enif_get_list_length(env, list, &nrowA))
        return enif_make_badarg(env);
    if(!enif_get_list_length(env, head, &ncolA))
        return enif_make_badarg(env);

    if (nrowA==0)
        return enif_make_badarg(env);

    double *a = enif_alloc(ncolA*nrowA * sizeof *a);
    ERL_NIF_TERM *A = enif_alloc(ncolA*nrowA * sizeof *A);

    double *pcoeffa;
    pcoeffa=&a[0];

    ERL_NIF_TERM *pcoeffA;
    pcoeffA=&A[0];

    while (enif_get_list_cell(env, list, &head, &tail)) {

        sublist=head;
        while (enif_get_list_cell(env, sublist, &subhead, &subtail)) {
            if (!enif_get_double(env, subhead, &val)) 
				return enif_make_badarg(env);

            *pcoeffa = val;
            pcoeffa++;

            *pcoeffA =enif_make_double(env,val);
            pcoeffA++;

            sublist = subtail;
        }
        list = tail;
    }
    *pcoeffa = '\0';
    *pcoeffA = '\0';

    list=argv[1];

    enif_get_list_cell(env, list, &head, &tail);
    if (!enif_get_list_length(env, list, &nrowB))
        return enif_make_badarg(env);
    if(!enif_get_list_length(env, head, &ncolB))
        return enif_make_badarg(env);

    if (nrowB==0)
        return enif_make_badarg(env);

    double *b = enif_alloc(ncolB*nrowB * sizeof *b);
    ERL_NIF_TERM *B = enif_alloc(ncolB*nrowB * sizeof *B);

    double *pcoeffb;
    pcoeffb=&b[0];

   ERL_NIF_TERM *pcoeffB;
    pcoeffB=&B[0];

    while (enif_get_list_cell(env, list, &head, &tail)) {

        sublist=head;
        while (enif_get_list_cell(env, sublist, &subhead, &subtail)) {
            if (!enif_get_double(env, subhead, &val))
				return enif_make_badarg(env);

            *pcoeffb = val;
            pcoeffb++;

            *pcoeffB =enif_make_double(env,val);
            pcoeffB++;

           sublist = subtail;
        }
        list = tail;
    }
    *pcoeffb = '\0';
    *pcoeffB = '\0';

//ERL_NIF_TERM ncolA_nif = enif_make_int(env, ncolA);
//ERL_NIF_TERM nrowA_nif = enif_make_int(env, nrowA);
//ERL_NIF_TERM ncolB_nif = enif_make_int(env, ncolB);
//ERL_NIF_TERM nrowB_nif = enif_make_int(env, nrowB);

//return enif_make_tuple2(env,enif_make_tuple3(env,enif_make_tuple2(env,enif_make_atom(env,"ncolA"),ncolA_nif),enif_make_tuple2(env,enif_make_atom(env,"nrowA"),nrowA_nif),enif_make_list_from_array(env, A, ncolA*nrowA)),enif_make_tuple3(env,enif_make_tuple2(env,enif_make_atom(env,"ncolB"),ncolB_nif),enif_make_tuple2(env,enif_make_atom(env,"nrowB"),nrowB_nif),enif_make_list_from_array(env, B, ncolB*nrowB)));

    int i,j;

    int nrowC=nrowA,ncolC=ncolB;
    double *c = enif_alloc(ncolB*nrowA * sizeof *c);
    for (i=0; i<nrowC; i++) {
    	for (j=0; j<ncolC; j++) {
		 	c[(i*ncolA)+j]=0;
		}
 	}

	int CblasRowMajor=101;
	int CblasNoTrans=111;
	cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
            nrowA,ncolB,ncolA,
            1.0,a,ncolA,
            b,ncolB,
            0.0,c,ncolC);

    ERL_NIF_TERM *C = enif_alloc(ncolB*nrowA * sizeof *C);
	ERL_NIF_TERM c_nif = enif_make_list(env, 0);

    for (i=nrowC-1; i>=0; i--) {
    	for (j=0; j<ncolC; j++) {
		 	C[j]=enif_make_double(env,c[(i*ncolC)+j]);
		}
		c_nif=enif_make_list_cell(env,enif_make_list_from_array(env, C, ncolC),c_nif);
 	}

    // and then clean them up too:
    enif_free(a);
    enif_free(b);
    enif_free(c);

    enif_free(A);
    enif_free(B);
    enif_free(C);

	return c_nif;


}

static ERL_NIF_TERM svd_nif(ErlNifEnv* env, int argc, const ERL_NIF_TERM argv[])
{

	unsigned int ncol=-1,nrow=-1;
	double val;

	ERL_NIF_TERM list, head, tail, sublist, subhead,subtail; 
	list=argv[0];
    enif_get_list_cell(env, list, &head, &tail);
	if (!enif_get_list_length(env, list, &ncol))
		return enif_make_badarg(env);
	if(!enif_get_list_length(env, head, &nrow))
		return enif_make_badarg(env);

	if (ncol==0)
		return enif_make_badarg(env);

    double *a = enif_alloc(ncol*nrow * sizeof *a);
    ERL_NIF_TERM *A = enif_alloc(ncol*nrow * sizeof *A);

	double *pcoeff;
	pcoeff=&a[0];

	ERL_NIF_TERM *pcoeffA;
	pcoeffA=&A[0];

    while (enif_get_list_cell(env, list, &head, &tail)) {
		sublist=head;
    	while (enif_get_list_cell(env, sublist, &subhead, &subtail)) {
			if (!enif_get_double(env, subhead, &val)) return enif_make_badarg(env);

			*pcoeff = val;  
			pcoeff++;  

		 	*pcoeffA =enif_make_double(env,val);
		 	pcoeffA++;

			sublist = subtail;   
    	}
		list = tail;   
	}
	*pcoeff = '\0';
	*pcoeffA = '\0';


//ERL_NIF_TERM ncol_nif = enif_make_int(env, ncol);
//ERL_NIF_TERM nrow_nif = enif_make_int(env, nrow);

//return enif_make_tuple3(env,enif_make_tuple2(env,enif_make_atom(env,"ncol"),ncol_nif),enif_make_tuple2(env,enif_make_atom(env,"nrow"),nrow_nif),enif_make_list_from_array(env, A, ncol*nrow));

  int m = nrow, n = ncol, lda = m, ldu = m, ldvt = n;
  /* Local arrays */
  /*double s[N], u[LDU*M], vt[LDVT*N];*/
//  double a[6*4] = {
//            7.52, -0.76,  5.13, -4.75,  1.33, -2.40,
//           -1.10,  0.62,  6.62,  8.52,  4.91, -6.77,
//           -7.95,  9.34, -5.66,  5.75, -5.49,  2.34,
//            1.08, -7.10,  0.87,  5.30, -3.52,  3.95
//  };

    // Setup a buffer to hold the singular values:
    int min_mn = m < n ? m : n;
    int max_mn= m > n ? m : n;
    double *s = enif_alloc(min_mn * sizeof *s);

    // Setup buffers to hold the matrices U and Vt:
    double *u = enif_alloc(m*m * sizeof *u);
    double *vt = enif_alloc(n*n * sizeof *vt);

    // Workspace and status variables:
    int *iwork = enif_alloc(8*min_mn);

    int lwork = 3*min_mn*min_mn+MAX(max_mn,4*min_mn*min_mn+4*min_mn);
    double *work = enif_alloc(lwork * sizeof *work);

    // Call dgesdd_ to do the actual computation:
    int info = 0;
    dgesdd_("A", &m, &n, a, &lda, s, u, &m, vt, &n, work, &lwork, iwork, &info);

    // Cleanup workspace:
    enif_free(work);
    enif_free(iwork);

    // do something useful with U, S, Vt ...
    int i,j;

	ERL_NIF_TERM *D = enif_alloc(min_mn * sizeof *D);
  	for (i=0; i<min_mn; ++i) D[i]=enif_make_double(env, s[i]);
	ERL_NIF_TERM d_nif = enif_make_list_from_array(env, D, min_mn);

	ERL_NIF_TERM *U = enif_alloc(m * sizeof *U);
	ERL_NIF_TERM u_nif = enif_make_list(env, 0);
    for (i=ldu-1; i>=0; i--) {
    	for (j=0; j<ldvt; j++) {
		 	U[j]=enif_make_double(env,u[i+j*ldu]);
		}
		u_nif=enif_make_list_cell(env,enif_make_list_from_array(env, U, ldvt),u_nif);
 	}

	ERL_NIF_TERM *V = enif_alloc(n * sizeof *V);
	ERL_NIF_TERM v_nif = enif_make_list(env, 0);
    for (i=ldvt-1; i>=0; i--) {
    	for (j=0; j<ldvt; j++) {
		 	V[j]=enif_make_double(env,vt[j+i*ldvt]);
		}
		v_nif=enif_make_list_cell(env,enif_make_list_from_array(env, V, ldvt),v_nif);
 	}

    // and then clean them up too:
    enif_free(s);
    enif_free(u);
    enif_free(vt);

	enif_free(D);
	enif_free(U);
	enif_free(V);

	// Note we switch the u and v lapack to match the R SVD convenvtion of left (v) and right (u)
    if (info>0) 
		return enif_make_tuple2(env,enif_make_atom(env,"error"),enif_make_string(env,  "Can not converge to answer", ERL_NIF_LATIN1));
	else if (info<0) 
  		return enif_make_badarg(env);
	else 
  		return enif_make_tuple3(env,enif_make_tuple2(env,enif_make_atom(env,"d"),d_nif),enif_make_tuple2(env,enif_make_atom(env,"u"),v_nif),enif_make_tuple2(env,enif_make_atom(env,"v"),u_nif));
}

static ErlNifFunc nif_funcs[] = {
    {"transpose", 1, transpose_nif},
    {"matmul", 2, matmul_nif},
    {"svd", 1, svd_nif}
};

ERL_NIF_INIT(linalg_lapack, nif_funcs, NULL, NULL, NULL, NULL)
