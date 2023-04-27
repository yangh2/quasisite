#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#define swap(a,b) {a=a+b;b=a-b;a=a-b;}

typedef struct point_type{
  int *site;
  struct point_type *prev;
  struct point_type *next;
  struct point_type *father;
  int idx;
} point_t;
typedef point_t *point_p_t;

int AllocPoint(point_t *point, int dim){
  point->site = (int *) malloc(sizeof(int)* dim);
  point->prev = NULL; point->next = NULL; point->father = NULL;
  point->idx = 0;
  return 0;
}
int IsSame(int dim, point_t *point1, point_t *point2){
  for (int i=0;i<dim;i++)
    if ( point1->site[i]!= point2->site[i])
      return 0;
  return 1;
}

double CalDistance(int dim, double *point, double *vector1, double *vector2){
  double a1=0, a2=0;
  for (int i=0;i<dim;i++) a1 += point[i]*vector1[i];
  for (int i=0;i<dim;i++) a2 += point[i]*vector2[i];
  double dis = 0;
  for (int i=0;i<dim;i++) dis += (point[i]-a1*vector1[i]-a2*vector2[i])*(point[i]-a1*vector1[i]-a2*vector2[i]);
  dis = sqrt(dis);
  return dis;
}

double CalMaxDistance(int dim, double *vector1, double *vector2){
  double distance = 0, tmp_distance;
  double *b;
  b = (double *) malloc(sizeof(double)*dim);
  for (int i=0;i<dim;i++) b[i] = 1;
  while (b[0] >= 0){
    tmp_distance = CalDistance(dim, b, vector1, vector2);
    distance = distance>tmp_distance?distance:tmp_distance;
    b[dim-1] -=1;
    for (int i=dim-1;i>0;i--)
      if ( b[i] < 0 ){
	b[i]=1;
	b[i-1]-=1;
      }
      else
	break;
  }
  return distance;
}

int Criteria1(int hdim, int pdim, point_t *point, point_t *center, double **basis, double *site, double radius){
  double *b; double d1=0,d2=0;
  b = (double *) malloc(sizeof(double)*hdim);
  for (int i=0;i<hdim;i++) b[i] = center->site[i]-site[i];
  for (int i=0;i<hdim;i++) d1+= b[i]*b[i]; d1=sqrt(d1);
  for (int i=0;i<hdim;i++) b[i] = point->site[i]-site[i];
  for (int i=0;i<hdim;i++) d2+= b[i]*b[i]; d2=sqrt(d2);
  double *a; a = (double *) malloc(sizeof(double)*hdim);
  for (int i=0;i<hdim;i++) a[i] = 0;
  for (int i=0;i<hdim;i++){
    for (int j=0;j<hdim;j++) a[i] += basis[i][j]*b[j];
  }
  double r1 = 0, dist = 0;
  for (int i=0;i<pdim;i++) r1 += a[i]*a[i]; r1 = sqrt(r1);
  for (int i=pdim;i<hdim;i++) dist += a[i]*a[i]; dist = sqrt(dist);
  free(b);
  //if ( (dist<2) && (r1<radius) && (d2>d1)) return 1;
  //if ((d2>d1)) return 1;
  if ( (r1<radius) && (d2>d1)) return 1;
  else return 0;
}

int Criteria2(int dim, point_t *point, double *vector1, double *vector2, double *site, double maxdis){
  double *b;
  b = (double *) malloc(sizeof(double)*dim);
  for (int i=0;i<dim;i++) b[i] = point->site[i]-site[i];
  double a1=0, a2=0;
  for (int i=0;i<dim;i++) a1 += b[i]*vector1[i];
  for (int i=0;i<dim;i++) a2 += b[i]*vector2[i];
  double dis = 0;
  for (int i=0;i<dim;i++) dis += (b[i]-a1*vector1[i]-a2*vector2[i])*(b[i]-a1*vector1[i]-a2*vector2[i]);
  dis = sqrt(dis);
  free(b);
  if  (dis<=maxdis) return 1;
  else return 0;
}
int Criteria3(int dim, point_t *point1, point_t *point2, double *vector1, double *vector2, double *site){
  double *b1, *b2;
  b1 = (double *) malloc(sizeof(double)*dim);
  b2 = (double *) malloc(sizeof(double)*dim);
  double a1=0, a2=0; double d1=0,d2=0;

  for (int i=0;i<dim;i++) b1[i] = point1->site[i]-site[i];
  for (int i=0;i<dim;i++) a1 += b1[i]*vector1[i];
  for (int i=0;i<dim;i++) a2 += b1[i]*vector2[i];
  for (int i=0;i<dim;i++) b1[i] = b1[i]-a1*vector1[i]-a2*vector2[i];
  for (int i=0;i<dim;i++) d1 += b1[i]*b1[i]; d1 = sqrt(d1);
  
  a1 = 0; a2= 0;
  for (int i=0;i<dim;i++) b2[i] = point2->site[i]-site[i];
  for (int i=0;i<dim;i++) a1 += b2[i]*vector1[i];
  for (int i=0;i<dim;i++) a2 += b2[i]*vector2[i];
  for (int i=0;i<dim;i++) b2[i] = b2[i]-a1*vector1[i]-a2*vector2[i];
  for (int i=0;i<dim;i++) d2 += b2[i]*b2[i]; d2 = sqrt(d2);

  a1 = 0;
  for (int i=0;i<dim;i++) a1 += b1[i]*b2[i];
  free(b1); free(b2);
  if ( fabs(a1/d1/d2 + 1)<0.00001 ) return 1; else return 0;
}
int Criteria5(int dim, point_t *point, double *vector1, double *vector2, double *site){
  double *b;
  b = (double *) malloc(sizeof(double)*dim);
  for (int i=0;i<dim;i++) b[i] = point->site[i]-site[i];
  double a1=0, a2=0;
  for (int i=0;i<dim;i++) a1 += b[i]*vector1[i];
  for (int i=0;i<dim;i++) a2 += b[i]*vector2[i];
  for (int i=0;i<dim;i++) b[i] -= (a1*vector1[i]+a2*vector2[i]);
  double dist =0 ; 
  for (int i=0;i<dim;i++) dist += b[i]*b[i]; dist = sqrt(dist);
  for (int i=0;i<dim;i++) b[i] = fabs(b[i])/dist;
  double dismax = 0;
  //for (int i=0;i<dim;i++) dismax += b[i]>0?b[i]:-b[i]; dismax /= 2;
  for (int i=0;i<dim;i++)
    for (int j=dim-1;j>i;j--)
      if ( b[j]>b[j-1] )
	swap(b[j], b[j-1]);
  dismax = b[0]+b[1];
  free(b);
  if ( (dist<dismax)) return 1;
  else return 0;
}
int Criteria6(int hdim, int pdim, point_t *point, double **basis, double *site){
  double *b;
  b = (double *) malloc(sizeof(double)*hdim);
  for (int i=0;i<hdim;i++) b[i] = point->site[i]-site[i];
  double *a; a = (double *) malloc(sizeof(double)*hdim);
  for (int i=0;i<hdim;i++) a[i] = 0;
  for (int i=0;i<hdim;i++){
    for (int j=0;j<hdim;j++) a[i] += basis[i][j]*b[j];
  }
  double fmax; fmax = 0;
  for (int i=pdim;i<hdim;i++) fmax = fmax>fabs(a[i])?fmax:fabs(a[i]);
  free(b); free(a);
  if ( fmax < 0.5 ) return 1;
  else return 0;
}

int Iteration(int hdim, int pdim, double **basis, double *b, int *idx, int oneloc){
  double **mat;
  mat = (double **) malloc(sizeof(double*)*(hdim - pdim));
  for (int i=0;i<hdim-pdim;i++) mat[i] = (double *) malloc(sizeof(double)*(hdim - pdim));

  double *a; a = (double *) malloc(sizeof(double)*hdim);
  double *vec; vec = (double *) malloc(sizeof(double)*hdim);
  double *translate; translate = (double *) malloc(sizeof(double)*hdim);
  double *translateb; translateb = (double *) malloc(sizeof(double)*hdim);

  gsl_matrix *m_inv = gsl_matrix_alloc(hdim-pdim, hdim-pdim);
  gsl_permutation * p = gsl_permutation_alloc (hdim-pdim);
  int sign;

  int count =0;
  gsl_matrix *m = gsl_matrix_alloc(hdim - pdim, hdim - pdim);
  gsl_vector *res = gsl_vector_alloc(hdim - pdim);
  gsl_vector *sol = gsl_vector_alloc(hdim - pdim);
  int flag = 1; int bcm;

  count = 0;
  for (int k=0;k<hdim;k++){
    if ( idx[k] == 0 ) continue;
    for (int i=0;i<hdim;i++) vec[i] = 0;
    vec[k]+=1 ;
    for (int i=0;i<hdim;i++) a[i] = 0;
    for (int i=0;i<hdim;i++){
      for (int j=0;j<hdim;j++) a[i] += basis[i][j]*vec[j];
    }
    for (int i=pdim;i<hdim;i++) mat[count][i-pdim] = a[i];
    count ++;
    vec[k]-= 1;
  }
  for (int i=0;i<hdim-pdim;i++)
    for (int j=0;j<hdim-pdim;j++)
      gsl_matrix_set(m, i, j, mat[i][j]);

  gsl_linalg_LU_decomp(m, p, &sign);
  if (fabs(gsl_linalg_LU_det(m,sign)) > 1e-10){
    gsl_linalg_LU_invert(m, p, m_inv);
    
    for (int bc=0;bc<(1<<(pdim));bc++){
      flag = 1;
      for (int i=0;i<hdim;i++) translate[i] = -0.5;
      bcm = bc;
      for (int k=0;k<hdim;k++){
	if ( idx[k] ==0 ) {
	  if (bcm % 2 == 1) translate[k] += 1;
	  bcm = bcm>>1;
	}
      }
      
      for (int i=0;i<hdim;i++) translateb[i]=0;
      for (int i=0;i<hdim;i++){
	for (int j=0;j<hdim;j++) translateb[i] += basis[i][j]*translate[j];
      }
      for (int i=0;i<hdim-pdim;i++) gsl_vector_set(res, i, b[i+pdim]-translateb[i+pdim]);
      
      gsl_blas_dgemv(CblasTrans, 1.0, m_inv, res, 0, sol);
      
      for (int i=0;i<hdim-pdim;i++) if ((gsl_vector_get(sol, i)>=1)||(gsl_vector_get(sol,i)<=0)) {flag = 0; break;}
      if (flag) break;
    }
  }
  else {flag=0;}
  
  free(a); free(vec);
  for (int i=0;i<hdim-pdim;i++) free(mat[i]);
  free(mat);
  gsl_matrix_free(m); gsl_vector_free(res);
  gsl_matrix_free(m_inv); gsl_permutation_free(p);
  if ( flag ) return 1;
  else {
    oneloc --;
    if ( oneloc < 0) return 0;
    for (int i=oneloc+1; i<hdim; i++){
      if ( idx[i] == 1 ) return 0;
      idx[i] = 1; idx[oneloc] = 0;
      flag = Iteration(hdim, pdim, basis, b, idx, oneloc);
      idx[i] = 0; idx[oneloc] = 1;
      if (flag) return flag;
    }
  }
  return 0;
}

int Criteria4(int hdim, int pdim, point_t *point, double **basis, double *site){
  double *b;
  b = (double *) malloc(sizeof(double)*hdim);
  for (int i=0;i<hdim;i++) b[i] = point->site[i]-site[i];
  double *a; a = (double *) malloc(sizeof(double)*hdim);
  for (int i=0;i<hdim;i++) a[i] = 0;
  for (int i=0;i<hdim;i++){
    for (int j=0;j<hdim;j++) a[i] += basis[i][j]*b[j];
  }
  int *idx;
  idx = (int *) malloc(sizeof(int)*hdim);
  for (int i=0;i<hdim-pdim;i++) idx[i] = 1; for ( int i=hdim-pdim;i<hdim;i++) idx[i] = 0;

  int flag = Iteration(hdim, pdim, basis, a, idx, hdim-pdim);
  //double fmax; fmax = 0;
  //for (int i=pdim;i<hdim;i++) fmax = fmax>fabs(a[i])?fmax:fabs(a[i]);
  free(b); free(a); free(idx);
  return flag;
}

int Normalize(int hdim, int pdim, double **basis){
  double **cart;
  cart = (double **) malloc(sizeof(double*)*hdim);
  for (int i=0;i<hdim;i++) cart[i] = (double *) malloc(sizeof(double)*hdim);
  for (int i=0;i<hdim;i++){
    for (int j=0;j<hdim;j++)
      cart[i][j] = 1;
    cart[i][i]=-1;
  }
  double dist = 0, prod = 0;
  for (int i=0;i<pdim;i++){
    dist = 0;
    for (int j=0;j<hdim;j++) dist += basis[i][j]*basis[i][j]; dist = sqrt(dist);
    for (int j=0;j<hdim;j++) basis[i][j] /= dist;
    for (int k=i+1;k<pdim;k++){
      prod = 0;
      for (int j=0;j<hdim;j++) prod += basis[i][j] * basis[k][j];
      for (int j=0;j<hdim;j++) basis[k][j] -= prod *basis[i][j];
    }
  }
  for (int i=0;i<pdim;i++){
    for (int k=0;k<hdim;k++){
      prod = 0;
      for (int j=0;j<hdim;j++) prod += basis[i][j] * cart[k][j];
      for (int j=0;j<hdim;j++) cart[k][j] -= prod *basis[i][j];
    }
  }
  for (int i=0;i<hdim;i++){
    dist = 0;
    for (int j=0;j<hdim;j++) dist += cart[i][j]*cart[i][j]; dist = sqrt(dist);
    if ( dist > 1e-10 ) for (int j=0;j<hdim;j++) cart[i][j] /= dist;
    for (int k=i+1;k<hdim;k++){
      prod = 0;
      for (int j=0;j<hdim;j++) prod += cart[i][j] * cart[k][j];
      for (int j=0;j<hdim;j++) cart[k][j] -= prod *cart[i][j];
    }
  }
  int idx = pdim;
  for (int i=0;i<hdim;i++){
    dist = 0;
    for (int j=0;j<hdim;j++) dist += cart[i][j]*cart[i][j]; dist = sqrt(dist);
    if ( dist > 1e-10 ){
      for (int j=0;j<hdim;j++) basis[idx][j] = cart[i][j];
      idx ++;
      if ( idx == hdim)  break;
    }
  }
  for (int i=0;i<hdim;i++){
    dist = 0;
    for (int j=0;j<hdim;j++) dist += basis[i][j]*basis[i][j]; dist = sqrt(dist);
    if ( dist > 1e-10 ) for (int j=0;j<hdim;j++) basis[i][j] /= dist;
    for (int k=i+1;k<hdim;k++){
      prod = 0;
      for (int j=0;j<hdim;j++) prod += basis[i][j] * basis[k][j];
      for (int j=0;j<hdim;j++) basis[k][j] -= prod *basis[i][j];
    }
  }
  for (int i=0;i<hdim;i++)
    for (int j=i+1;j<hdim;j++) swap(basis[i][j], basis[j][i]);

  gsl_matrix *m = gsl_matrix_alloc(hdim, hdim);
  gsl_matrix *m_inv = gsl_matrix_alloc(hdim, hdim);
  for (int i=0;i<hdim;i++)
    for (int j=0;j<hdim;j++)
      gsl_matrix_set(m, i, j, basis[i][j]);
  gsl_permutation * p = gsl_permutation_alloc (hdim);
  int sign;
  
  gsl_linalg_LU_decomp(m, p, &sign);
  gsl_linalg_LU_invert(m, p, m_inv);
  for (int i=0;i<hdim;i++)
    for (int j=0;j<hdim;j++)
      basis[i][j] = gsl_matrix_get(m_inv, i, j);

  gsl_matrix_free(m); gsl_matrix_free(m_inv);
  gsl_permutation_free(p);
  return 0;
}

int main(){
  int hdim, pdim;
  double **bvectors;
  FILE *fid = fopen("BVECTORS", "r");
  fscanf(fid, "%d %d", &hdim, &pdim);
  bvectors = (double **) malloc(sizeof(double *)*hdim);
  for (int i=0;i<hdim;i++) bvectors[i] = (double *) malloc(sizeof(double)*pdim);
  for (int i=0;i<hdim;i++)
    for (int j=0;j<pdim;j++)
      fscanf(fid, "%lf", &(bvectors[i][j]));
  fclose(fid);

  int hdim1;
  fid = fopen("PLANE", "r");
  fscanf(fid, "%d", &hdim1);
  if (hdim1 != hdim) {printf("inconsistent dimensions in BVECTORS and PLANE!\n"); return (0);}
  double **basis;
  basis = (double **) malloc(sizeof( double *)*hdim);
  for (int i=0;i<hdim;i++) basis[i] = (double *) malloc(sizeof(double)*hdim);
  for (int i=0;i<pdim;i++)
    for ( int j=0;j<hdim;j++)
      fscanf(fid, "%lf", &(basis[i][j]));
  /* for (int i=0;i<hdim;i++){ */
  /*   for (int j=0;j<hdim;j++) */
  /*     printf("%lf ", basis[i][j]); */
  /*   printf("\n"); */
  /* } */
  Normalize(hdim, pdim, basis);
  for (int i=0;i<hdim;i++){
    for (int j=0;j<hdim;j++)
      printf("%lf ", basis[i][j]);
    printf("\n");
  }
  
  double *site;
  site = (double *) malloc(sizeof(double)*hdim);
  for (int i=0;i<hdim;i++) fscanf(fid, "%lf", &(site[i]));
  double radius;
  fscanf(fid, "%lf", &radius);
  fclose(fid);

  point_t *head=(point_t *) malloc(sizeof(point_t));  AllocPoint(head, hdim);
  point_t *point=(point_t *) malloc(sizeof(point_t));  AllocPoint(point, hdim);
  for (int i=0;i<hdim;i++) point->site[i] = round(site[i]);
  head->next = point; point->prev = head;

  int flag= 0; point_t *ptr; int free_flag; int min_flag;
  fid = fopen("points.dat", "w");
  point_t *tail=point, *center; int sqr_flag1, sqr_flag2;
  while (head->next != NULL){
    center = head->next; head->next=center->next;
    if ( center->next != NULL) center->next->prev = head;
    if ( tail == center ) tail = head;
    for (int iter=0;iter<hdim;iter++){
      min_flag = 1;
      point=(point_t *) malloc(sizeof(point_t));  AllocPoint(point, hdim); free_flag=1;
      for (int j=0;j<hdim;j++) point->site[j] = center->site[j];
      point->site[iter] += 1; point->idx = iter;
      if ( Criteria1(hdim, pdim, point, center, basis, site, radius) ){
	if ( Criteria4(hdim, pdim,  point, basis, site) ){
	  for (int k=0;k<hdim;k++){
	    if ( k== iter ) continue;
	    sqr_flag1 = 0; sqr_flag2 = 0;
	    point->site[k] += 1;
	    if ( Criteria4(hdim, pdim, point, basis,  site) ) sqr_flag1 = 1;
	    point->site[iter] -= 1;
	    if ( Criteria4(hdim, pdim, point, basis,  site) ) sqr_flag2 = 1;
	    point->site[iter] += 1; point->site[k] -= 1;

	    if ( sqr_flag1 && sqr_flag2 ){
	      min_flag = 0;
	      for (int i=0;i<hdim;i++) fprintf(fid, "%d ", center->site[i]);
	      for (int i=0;i<hdim;i++) fprintf(fid, "%d ", point->site[i]);
	      fprintf(fid,"\n");
	      flag = 0; ptr = tail;
	      while (ptr->prev != NULL) {
		if ( IsSame(hdim, point, ptr) ) { flag=1; break;}
		ptr=ptr->prev;
	      }
	      if ( !flag ){
		point->next= tail->next; tail->next = point; point->prev = tail;
		point->father = center; tail = tail->next;
		free_flag = 0;
	      }
	      break;
	    }
	  }
	  if (min_flag == 1){
	    for (int k=0;k<hdim;k++){
	      if ( k== iter ) continue;
	      sqr_flag1 = 0; sqr_flag2 = 0;
	      point->site[k] -= 1;
	      if ( Criteria4(hdim, pdim, point, basis, site) ) sqr_flag1 = 1;
	      point->site[iter] -= 1;
	      if ( Criteria4(hdim, pdim, point, basis, site) ) sqr_flag2 = 1;
	      point->site[iter] += 1; point->site[k] += 1;

	      if ( sqr_flag1 && sqr_flag2 ){
		for (int i=0;i<hdim;i++) fprintf(fid, "%d ", center->site[i]);
		for (int i=0;i<hdim;i++) fprintf(fid, "%d ", point->site[i]);
		fprintf(fid,"\n");
		flag = 0; ptr = tail;
		while (ptr->prev != NULL) {
		  if ( IsSame(hdim, point, ptr) ) { flag=1; break;}
		  ptr=ptr->prev;
		}
		if ( !flag ){
		  point->next= tail->next; tail->next = point; point->prev = tail;
		  point->father = center; tail = tail->next;
		  free_flag = 0;
		}
		break;
	      }
	    }
	  }
	}
      }
      if (free_flag ) free(point);

      min_flag = 1;
      point=(point_t *) malloc(sizeof(point_t));  AllocPoint(point, hdim); free_flag = 1;
      for (int j=0;j<hdim;j++) point->site[j] = center->site[j];
      point->site[iter] -= 1; point->idx = iter;
      if ( Criteria1(hdim, pdim, point, center, basis,site, radius) ){
	if ( Criteria4(hdim, pdim, point, basis, site) ) {
	  for (int k=0;k<hdim;k++){
	    if ( k== iter ) continue;
	    sqr_flag1 = 0; sqr_flag2 = 0;
	    point->site[k] -= 1;
	    if ( Criteria4(hdim, pdim, point, basis, site) ) sqr_flag1 = 1;
	    point->site[iter] += 1;
	    if ( Criteria4(hdim, pdim, point, basis, site) ) sqr_flag2 = 1;
	    point->site[iter] -= 1; point->site[k] += 1;

	    if ( sqr_flag1 && sqr_flag2 ){
	      min_flag = 0;
	      for (int i=0;i<hdim;i++) fprintf(fid, "%d ", center->site[i]);
	      for (int i=0;i<hdim;i++) fprintf(fid, "%d ", point->site[i]);
	      fprintf(fid,"\n");
	      flag = 0; ptr = tail;
	      while (ptr->prev != NULL) {
		if ( IsSame(hdim, point, ptr) ) { flag=1; break;}
		ptr=ptr->prev;
	      }
	      if ( !flag ){
		point->next= tail->next; tail->next = point; point->prev = tail;
		point->father = center; tail = tail->next;
		free_flag = 0;
	      }
	      break;
	    }
	  }
	  if (min_flag == 1){
	    for (int k=0;k<hdim;k++){
	      if ( k== iter ) continue;
	      sqr_flag1 = 0; sqr_flag2 = 0;
	      point->site[k] += 1;
	      if ( Criteria4(hdim, pdim, point, basis, site) ) sqr_flag1 = 1;
	      point->site[iter] += 1;
	      if ( Criteria4(hdim, pdim, point, basis, site) ) sqr_flag2 = 1;
	      point->site[iter] -= 1; point->site[k] -= 1;

	      if ( sqr_flag1 && sqr_flag2 ){
		for (int i=0;i<hdim;i++) fprintf(fid, "%d ", center->site[i]);
		for (int i=0;i<hdim;i++) fprintf(fid, "%d ", point->site[i]);
		fprintf(fid,"\n");

		flag = 0; ptr = tail;
		while (ptr->prev != NULL) {
		  if ( IsSame(hdim, point, ptr) ) { flag=1; break;}
		  ptr=ptr->prev;
		}
		if ( !flag ){
		  point->next= tail->next; tail->next = point; point->prev = tail;
		  point->father = center; tail = tail->next;
		  free_flag = 0;
		}
		break;
	      }
	    }
	  }
	}
      }
      if (free_flag ) free(point);
    }
    free(center);
  }
  fclose(fid);
  for (int i=0;i<hdim;i++) free(bvectors[i]);
  free(bvectors);
  for (int i=0;i<hdim;i++) free(basis[i]);
  free(basis);
}

