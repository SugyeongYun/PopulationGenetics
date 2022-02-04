
#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#define num_ind (50)
#define num_gen (1000)
#define mut_rate (0.0005)
#define rounds (300)
#define seq_len (1)
#define K (50)
#define g (20)


long idum_num = 65248555; //random number
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


double get_het(double count_A, double count_a) {
	double freq_A = count_A / (count_A + count_a);
	double freq_a = count_a / (count_A + count_a);
	return 1 - (freq_A * freq_A) - (freq_a * freq_a);
}
int mutation(int allele) {

	int prob = 1 / mut_rate;

	if (rand() % prob == 0) {
		return abs(allele - 1);
	}

	else {
		return allele;
	}
}



int main(void) {

	FILE* fp = fopen("poi+mut+seq.txt", "w");

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


		double** het = malloc(sizeof(double) * num_gen);
		for (int i = 0; i < num_gen; ++i) {
			het[i] = malloc(sizeof(double) * seq_len);
		}

		int N = num_ind;
		
		
		int heti = 0;
		for (int j = 0; j < seq_len; ++j) {
				het[heti][j] = get_het(num_ind, 0);
				fprintf(fp, "%lf,", het[heti][j]);
		}
		

		for (int j = 1; j <= num_gen - 1; ++j) {
			int** new_gen = malloc(sizeof(int) * num_ind * 2);
			for (int i = 0; i < num_ind * 2; ++i) {
				new_gen[i] = malloc(sizeof(int) * seq_len);
			}

			double* count_a = malloc(sizeof(double) * seq_len);
			for (int l = 0; l<seq_len; ++l) {
				count_a[l] = 0;
			}
			int idx = 0;
			for (int k = 0; k <= N - 1; ++k) {
				double lambda = (double) (1 + g) / (double) (1 + g * N / K);
				int offspring = poidev(lambda);
				for (int x = idx; x < idx + offspring; ++x) {
					for (int l = 0; l <= seq_len - 1; ++l) {
						int new_allele = ori_gen[k][l];
						new_allele = mutation(new_allele);
						if (new_allele == 1) {
							++count_a[l];
						}
						new_gen[x][l] = new_allele;
					}
				}
				idx = idx + offspring;
			}
			N = idx;
			
			for (int l = 0; l <= seq_len - 1; ++l) {
				het[j][l] = get_het(N - count_a[l], count_a[l]);
				fprintf(fp, "%lf,", het[j][l]);
			}

			for (int i = 0; i < num_ind * 2; ++i) {
				for (int j = 0; j < seq_len; ++j) {
					ori_gen[i][j] = new_gen[i][j];
				}
			}

			for (int i = 0; i < num_ind * 2; ++i) {
				free(new_gen[i]);
			}

			free(new_gen);
			free(count_a);


		}
		for (int i = 0; i < num_ind * 2; ++i) {
			free(ori_gen[i]);
		}
		
		free(ori_gen);
		free(het);
	}

	fclose(fp);


	return 0;
}
