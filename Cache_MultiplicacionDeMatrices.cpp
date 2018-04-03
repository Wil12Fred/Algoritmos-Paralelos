#include <iostream>
#include <ctime> 

using namespace std;
const int n=1002,m=1002,l=1002;

unsigned t0, t1;

void matrix_mult_wiki_block(const float*A , const float* B, float* C,
							const int N, const int M, const int K,int b) {
	const int block_size = 32/sizeof(float);
	for(int i=0; i<N; i++) {
		for(int j=0; j<K; j++) {
			C[K*i + j] = 0;
		}
	}
	t0=clock();
	for (int i0 = 0; i0 < N; i0 += block_size) {
		int imax = i0 + block_size > N ? N : i0 + block_size;
		int mi,ki,kij,sj,jmax,kmax,j0,k0,j1;
		for (int j0 = 0; j0 < M; j0 += block_size) {
			int jmax = j0 + block_size > M ? M : j0 + block_size;
			
			for (int k0 = 0; k0 < K; k0 += block_size) {
				 kmax = k0 + block_size > K ? K : k0 + block_size;
				
				for (int j1 = j0; j1 < jmax; ++j1) {
					 sj = M * j1;
					
					for (int i1 = i0; i1 < imax; ++i1) {
						 mi = M * i1;
						 ki = K * i1;
						 kij = ki + j1;
						
						for (int k1 = k0; k1 < kmax; ++k1) {
							C[kij] += A[mi + k1] * B[sj + k1];
						}
					}
				}
			}
		}
	}
	t1=clock();
	double time = (double(t1-t0)/CLOCKS_PER_SEC);
	cout << "Execution Time: " << time << endl;
}

float M1[n][m], M2[m][l],M1I[m][n],M2I[l][m],Ans[n][l];

int main() {
	
	t0=clock();
	for (int i=0;i<n;i++){
		for (int j=0;j<m;j++){
			for (int k=0;k<l;k++){
				Ans[i][k]+=M1[i][j]*M2[j][k];
			}
		}
	}
	t1=clock();
	double time = (double(t1-t0)/CLOCKS_PER_SEC);
	cout << "Execution Time: " << time << endl;
	
	/*t0=clock();
	for (int i=0;i<n;i++){
		for (int j=0;j<m;j++){
			for (int k=0;k<l;k++){
				Ans[i][k]+=M1[i][j]*M2I[k][j];
			}
		}
	}
	t1=clock();
	time = (double(t1-t0)/CLOCKS_PER_SEC);
	cout << "Execution Time: " << time << endl;*/
	matrix_mult_wiki_block(&M1[0][0],&M2[0][0],&Ans[0][0],n,m,l,6);
	matrix_mult_wiki_block(&M1[0][0],&M2[0][0],&Ans[0][0],n,m,l,2);
	matrix_mult_wiki_block(&M1[0][0],&M2[0][0],&Ans[0][0],n,m,l,3);
	matrix_mult_wiki_block(&M1[0][0],&M2[0][0],&Ans[0][0],n,m,l,4);
	matrix_mult_wiki_block(&M1[0][0],&M2[0][0],&Ans[0][0],n,m,l,5);
	matrix_mult_wiki_block(&M1[0][0],&M2[0][0],&Ans[0][0],n/2,m/2,l/2,6);
	matrix_mult_wiki_block(&M1[0][0],&M2[0][0],&Ans[0][0],n/2,m/2,l/2,2);
	matrix_mult_wiki_block(&M1[0][0],&M2[0][0],&Ans[0][0],n/2,m/2,l/2,3);
	matrix_mult_wiki_block(&M1[0][0],&M2[0][0],&Ans[0][0],n/2,m/2,l/2,4);
	matrix_mult_wiki_block(&M1[0][0],&M2[0][0],&Ans[0][0],n/2,m/2,l/2,5);
	
	return 0;
}

