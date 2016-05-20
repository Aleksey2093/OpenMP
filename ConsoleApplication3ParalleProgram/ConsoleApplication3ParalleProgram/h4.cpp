#include <stdio.h>
#include <chan.h>

double getKoef(int i)
{
	if ((i + 1) % 2)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

double getFactorial(int f)
{
	int i;
	double fact = 0;
	for (i = 1; i < f + 1; i++)
	{
		fact += i;
	}
	return fact;
}

double getStepen(double x, int st)
{
	int i;
	for (i = 0; i < st; i++)
	{
		x = x * x;
	}
	return x;
}


main(argc, argv, envp, in_ports, ins, out_ports, outs)
int argc, ins, outs;
char *argv[], *envp[];
CHAN *in_ports[], *out_ports[];
{
	double x, koef, stepen, factorial;
	int tsum, sum = 0, chunk, size, i, j;

	chan_in_word(&chunk, in_ports[0]);
	chan_in_word(&x, in_ports[0]);

	for (i = 3*chunk, j = 0; j < chunk; i++, j++)
	{
		koef = getKoef(i + 1);
		stepen = getStepen(x, i + 1);
		factorial = getFactorial(i + 1);
		sum += koef*stepen*factorial;
	}

	chan_out_word(sum, out_ports[0]);
}