#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#define num_ind (100)
#define num_gen (5800)
#define rounds (1000)
#define seq_len (990)
#define num_bn (55)
#define K (100)
#define s (0.05)
#define mut_rate (0.000084)

long idum_num = 65248; //random number
long* idum = &idum_num;

#define PI 3.141592654 //from poidev()
#define IA 16807	   //from ran1()
#define IM 2147483647
#define AM (1.0 / IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1 + (IM - 1) / NTAB)
#define EPS 1.2e-7
#define RNMX (1.0 - EPS)


double gammln(double xx)
{ //from poidev()
	double x, y, tmp, ser;
	static double cof[6] = { 76.18009172947146, -86.50532032941677, 24.01409824083091, -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5 };
	int j;

	y = x = xx;
	tmp = x + 5.5;
	tmp -= (x + 0.5) * log(tmp);
	ser = 1.000000000190015;
	for (j = 0; j <= 5; j++)
	{
		ser += cof[j] / ++y;
	}

	return -tmp + log(2.5066282746310005 * ser / x);
}


double ran1()
{ //seed number�� �ѹ� �־��ָ�, ���� ran1() �Լ��� ȣ���� ������ ��� �ٸ� random number ��ȯ (uniform distn.�� ������ random number generator) -> ���� ���� seed number�� �־��ָ� ���� ������ ����
	int j;
	long k;
	static long iy = 0;
	static long iv[NTAB];
	double temp;

	if (*idum <= 0 || !iy)
	{
		if (-(*idum) < 1)
		{
			*idum = 1;
		}
		else
		{
			*idum = -(*idum);
		}

		for (j = NTAB + 7; j >= 0; j--)
		{
			k = (*idum) / IQ;
			*idum = IA * (*idum - k * IQ) - IR * k;
			if (*idum < 0)
			{
				*idum += IM;
			}
			if (j < NTAB)
			{
				iv[j] = *idum;
			}
		}
		iy = iv[0];
	}
	k = (*idum) / IQ;
	*idum = IA * (*idum - k * IQ) - IR * k;
	if (*idum < 0)
	{
		*idum += IM;
	}
	j = iy / NDIV;
	iy = iv[j];
	iv[j] = *idum;
	if ((temp = AM * iy) > RNMX)
	{
		return RNMX;
	}
	else
	{
		return temp;
	}
} //��ȯ�� ����: 0 < ran1(idum) < 1, 0 <= (int)(ran1(idum)*num) < num


double poidev(double xm)
{ //xm: mean, ran1(idum): source of uniform random deviates
	//double gammln(double xx);
	//double ran1(long *idum);
	static double sq, alxm, Mu, oldm = (-1.0);
	double em, t, y;

	if (xm < 12.0)
	{
		if (xm != oldm)
		{
			oldm = xm;
			Mu = exp(-xm);
		}
		em = -1;
		t = 1.0;
		do
		{
			++em;
			t *= ran1(idum);
		} while (t > Mu);
	}
	else
	{
		if (xm != oldm)
		{
			oldm = xm;
			sq = sqrt(2.0 * xm);
			alxm = log(xm);
			Mu = xm * alxm - gammln(xm + 1.0);
		}
		do
		{
			do
			{
				y = tan(PI * ran1(idum));
				em = sq * y + xm;
			} while (em < 0.0);
			em = floor(em);
			t = 0.9 * (1.0 + y * y) * exp(em * alxm - gammln(em + 1.0) - Mu);
		} while (ran1(idum) > t);
	}
	return em;
}


int largest(int arr[], int n)
{
    int i;
    
    // Initialize maximum element
    int max = arr[0];
 
    // Traverse array elements from second and
    // compare every element with current max 
    for (i = 1; i < n; i++)
        if (arr[i] > max)
            max = arr[i];
 
    return max;
}

int sum(int arr[], int n) {
    int sum_arr = 0;
    for (int i = 0; i < n; ++i) {
        sum_arr += arr[i];
    }

    return sum_arr;
}

int mutation(int allele) {

	if (ran1() < mut_rate) {
		return abs(allele - 1);
	}

	else {
		return allele;
	}
}

int oneway_mut(int allele) {

	if (ran1() < mut_rate) {
		return 1;
	}

	else {
		return allele;
	}
}


int main(void) {

	double g = 0.2;

	FILE* fp1 = fopen("g02_beneficial.txt", "w");
	FILE* fp2 = fopen("g02_neutral.txt", "w");
	FILE* fp3 = fopen("g02_N.txt", "w");


	for (int i = 0; i < rounds; ++i) {

		int substitution_b = 0;
		int substitution_n = 0;

		int** ori_gen = (int **)malloc(sizeof(int *) * num_ind * 10);
		for (int i = 0; i < num_ind * 10; ++i) {
			ori_gen[i] = (int *) malloc(sizeof(int) * seq_len);
		}

		for (int i = 0; i < num_ind * 10; ++i) {
			for (int j = 0; j < seq_len; ++j) {
				ori_gen[i][j] = 0;
			}
		}

		double** freq_A = (double **)malloc(sizeof(double *) * num_gen);
		for (int i = 0; i < num_gen; ++i) {
			freq_A[i] = (double *) malloc(sizeof(double) * seq_len);
		}

		double** freq_a = (double **)malloc(sizeof(double *) * num_gen);
		for (int i = 0; i < num_gen; ++i) {
			freq_a[i] = (double *) malloc(sizeof(double) * seq_len);
		}


		int* ori_bn = (int *) malloc(sizeof(int) * num_ind * 10);

		for (int i = 0; i < num_ind * 10; ++i) {
			ori_bn[i] = 0;
		}


		int N = num_ind;

		for (int j = 1; j < num_gen; ++j) {
			int** new_gen = (int **) malloc(sizeof(int *) * num_ind * 10);
			for (int i = 0; i < num_ind * 10; ++i) {
				new_gen[i] = (int *) malloc(sizeof(int) * seq_len);
			}

            int* new_bn = (int *) malloc(sizeof(int) * num_ind * 10);
			for (int i = 0; i < num_ind * 10; ++i) {
				new_bn[i] = 0;
			}

            int bnmax = largest(ori_bn, num_ind * 10);
            double wmax = pow(1+s, bnmax);


			int* count_A = (int *) malloc(sizeof(int) * seq_len);
			for (int l = 0; l < seq_len; ++l) {
				count_A[l] = 0;
			}

			int* count_a = (int *) malloc(sizeof(int) * seq_len);
			for (int l = 0; l < seq_len; ++l) {
				count_a[l] = 0;
			}


			int idx = 0;
			for (int k = 0; k < N; ++k) {
				double lambda = (double) ((1 + g) * pow(1+s, ori_bn[k])) / (double) (wmax * (1 + g * N / K));
				int offspring = poidev(lambda);
				for (int x = idx; x < idx + offspring; ++x) {
					for (int l = 0; l < num_bn; ++l) {
						int new_allele = ori_gen[k][l];
                        new_allele = oneway_mut(new_allele);
						if (new_allele == 1) {
							++new_bn[x];
						}
						new_gen[x][l] = new_allele;
					}
					for (int l = num_bn; l < seq_len; ++l) {
						int new_allele = ori_gen[k][l];
                        new_allele = mutation(new_allele);
						if (new_allele == 0) {
							++count_A[l];
						}
						else if (new_allele == 1) {
							++count_a[l];
						}
						new_gen[x][l] = new_allele;
					}

				}
				idx = idx + offspring;
			}
			N = idx;

			for (int l = num_bn; l < seq_len; ++l) {
				freq_A[j][l] = (double) count_A[l] / (double) N;
				freq_a[j][l] = (double) count_a[l] / (double) N;

			}

            for (int i = 0; i < num_bn; ++i) {
                int sum_vertical = 0;
                for (int j = 0; j < N; ++j) {
                    sum_vertical += new_gen[j][i];
                }
                if (sum_vertical == N && N > 0) {
                    ++substitution_b;
                    for (int j = 0; j < N; ++j) {
                        new_gen[j][i] = 0;
                    }
                }

            }


			if (j > 3999) {
				fprintf(fp1, "%d,", substitution_b);
				fprintf(fp3, "%d,", N);
			}

			for (int i = 0; i < num_ind * 10; ++i) {
				for (int j = 0; j < seq_len; ++j) {
					ori_gen[i][j] = new_gen[i][j];
				}
			}

			for (int i = 0; i < num_ind * 10; ++i) {
				ori_bn[i] = new_bn[i];
			}

			free(count_A);
			free(count_a);

			for (int i = 0; i < num_ind * 10; ++i) {
				free(new_gen[i]);
			}

			free(new_gen);
			free(new_bn);

		}

		
		int derived[seq_len];
		
		for (int i = num_bn; i < seq_len; ++i) {
			if (freq_A[i] >= freq_a[i]) {
				derived[i] = 1;
			}
			else if (freq_A[i] < freq_a[i]) {
				derived[i] = 0;
			}
		}

		for (int j = 0; j < num_gen; ++j) {
			for (int i = num_bn; i < seq_len; ++i) {
				if (derived[i] == 1) {
					if (freq_a[j][i] > 0.99) {
						++substitution_n;
						derived[i] = 0;
					}
				}
				else if (derived[i] == 0) {
					if (freq_A[j][i] > 0.99) {
						++substitution_n;
						derived[i] = 1;
					}
				}
			}
			if (j > 3999) {
				fprintf(fp2, "%d,", substitution_n);
			}
		}

		for (int i = 0; i < num_ind * 2; ++i) {
			free(ori_gen[i]);
		}

		free(ori_gen);
		free(ori_bn);

		for (int i = 0; i < num_gen; ++i) {
			free(freq_A[i]);
		}

		free(freq_A);

		for (int i = 0; i < num_gen; ++i) {
			free(freq_a[i]);
		}

		free(freq_a);

	}

	fclose(fp1);
	fclose(fp2);
	fclose(fp3);
	


	return 0;
}
