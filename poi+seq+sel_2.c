
#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#define num_ind (200)
#define rounds (300)
#define seq_len (5)
#define num_bn (1)
#define K (200)
#define g (20)
#define s (0.05)

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


int main(void) {

	FILE* fp = fopen("poi+seq+sel(while).txt", "w");

	for (int i = 0; i < rounds; ++i) {

		int** ori_gen = malloc(sizeof(int) * num_ind * 2);
		for (int i = 0; i < num_ind * 2; ++i) {
			ori_gen[i] = malloc(sizeof(int) * seq_len);
		}

		for (int i = 0; i < num_ind * 2; ++i) {
			for (int j = 0; j < seq_len; ++j) {
				ori_gen[i][j] = 0;
			}
		}

		int* ori_bn = malloc(sizeof(int) * num_ind * 2);

		for (int i = 0; i < num_ind * 2; ++i) {
			ori_bn[i] = 0;
		}

		for (int i = 0; i < num_bn; ++i) {
			int x = rand() % num_ind;
			ori_gen[x][i] = 1;
			++ori_bn[x];
		}


		int N = num_ind;
		
		int sum_bn = sum(ori_bn, N);

		while (1) {
			int** new_gen = malloc(sizeof(int) * num_ind * 2);
			for (int i = 0; i < num_ind * 2; ++i) {
				new_gen[i] = malloc(sizeof(int) * seq_len);
			}
			int* new_bn = malloc(sizeof(int) * num_ind * 2);
			for (int i = 0; i < num_ind * 2; ++i) {
				new_bn[i] = 0;
			}

			int bnmax = largest(ori_bn, num_ind * 2);
			int wmax = pow(1+s, bnmax);

			int idx = 0;
			for (int k = 0; k <= N - 1; ++k) {

				double lambda = (double) ((1 + g) * pow(1+s, ori_bn[k])) / (double) (wmax * (1 + g * N / K));
				int offspring = poidev(lambda);
				for (int x = idx; x < idx + offspring; ++x) {
					for (int l = 0; l <= seq_len - 1; ++l) {
						int new_allele = ori_gen[k][l];
						if (new_allele == 1) {
						}
						new_gen[x][l] = new_allele;
					}
					for (int i = 0; i < num_bn; ++i) {
						if (new_gen[x][i] == 1) {
							++new_bn[x];
						}
					}

				}
				idx = idx + offspring;
			}
			N = idx;
			
			for (int i = 0; i < num_ind * 2; ++i) {
				for (int j = 0; j < seq_len; ++j) {
					ori_gen[i][j] = new_gen[i][j];
				}
			}

			for (int i = 0; i < num_ind * 2; ++i) {
				ori_bn[i] = new_bn[i];
			}

            sum_bn = sum(ori_bn, N);


			for (int i = 0; i < num_ind * 2; ++i) {
				free(new_gen[i]);
			}

			free(new_gen);
			free(new_bn);
			
			if (sum_bn == 0) {
				break;
			} else if (sum_bn == N) {
				break;
			}


		}

		fprintf(fp, "%d,", sum_bn);
		for (int i = 0; i < num_ind * 2; ++i) {
			free(ori_gen[i]);
		}

		
		free(ori_gen);
		free(ori_bn);
	}

	fclose(fp);


	return 0;
}
