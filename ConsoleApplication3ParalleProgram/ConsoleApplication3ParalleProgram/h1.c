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
	int size, tsum, sum = 0, chunk, i;

	printf("Vvedite N: ");
	scanf("%d",&size);
	printf("Vvedite X: ");
	scanf("%d",&x);
	chunk = (size) / 4;

	chan_out_word(chunk, out_ports[0]);
	chan_out_word(chunk, out_ports[1]);
	chan_out_word(chunk, out_ports[2]);

	sum = x;
	for (i = 1; i < chunk; i++)
	{
		koef = getKoef(i + 1);
		stepen = getStepen(x, i + 1);
		factorial = getFactorial(i + 1);
		sum += koef*stepen*factorial;
	}

	for (i = 0; i < 3; i++)
	{
		chan_in_word(&tsum, in_ports[i]);
		sum += tsum;
	}
	printf("Sum = %f", x);
}