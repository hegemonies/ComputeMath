#include <stdio.h>
#include <limits.h>

void swap(float *arr, int orig, int dest, int k)
{
	float tmp;
	for (int i = 0; i < k; i++) {
		tmp = arr[orig * k + i];
		arr[orig * k + i] = arr[dest * k + i];
		arr[dest * k + i] = tmp;
	}
}

int main()
{	
	float a[][4] = {
		{0, -1, 1, -1},
		{2, 0, 2, 3},
		{3, -1, 1, 3}
	};

	printf("\nOrigin:\n");
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 4; j++) {
			printf("%.2f ", a[i][j]);
		}
		printf("\n");
	}

	int n = 3;

	for (int i = 0; i < n - 1; i++) {
		printf("\n");
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 4; j++) {
				printf("%.2f ", a[i][j]);
			}
			printf("\n");
		}

		for (int j = i + 1; j < n; j++) {
			if (a[i][i] == 0) {
				int tmp = i;
				while (a[tmp][i] == 0) {
					tmp++;
				}
				swap(a, i, tmp, 4);
				printf("\nSWAP:\n");
				for (int i = 0; i < 3; i++) {
					for (int j = 0; j < 4; j++) {
						printf("%.2f ", a[i][j]);
					}
					printf("\n");
				}
			}
			// float max = INT_MIN;
			// int k;
			// for (int g = i; g < n - 1; g++) {
			// 	if (a[g][i] > max) {
			// 		max = a[g][i];
			// 		k = g;
			// 	}
			// }
			// if (i + 1 > k) {
			// 	printf("SWAP\n");
			// 	swap(a, i, k, n + 1);
			// }

			float mult = - (a[j][i] / a[i][i]);
			// printf("%.2f * %.2f\n", a[j][i], a[i][i]);
			// printf("mult = %.2f\n", mult);

			for (int k = i; k < n + 1; k++) {
				a[j][k] += a[i][k] * mult;
			}
		}
	}

	printf("\n");
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 4; j++) {
			printf("%.2f ", a[i][j]);
		}
		printf("\n");
	}

	// float z = a[2][2] / a[2][3];
	// float y = - (z * a[1][2] - a[1][3]) / a[1][1];
	// float x = - ((z * a[0][2])  + (y * a[0][1]) - a[0][3]) / a[0][0];

	float arr[3] = { 0 };

	int count = n - 1;

	for (int i = n - 1; i >= 0; i--) {
		float tmp = 0;
		int j;
		for (j = n - 1; j >= i; j--) {
			if (arr[j] != 0) {
				tmp += a[i][j] * arr[j];
				//printf("tmp = %.2f\n", tmp);
			}
		}
		arr[count] = (a[i][n] - tmp) / a[i][j + 1];
		count--;
	}

	for (int i = 0; i < n; i++) {
		printf("arr[%d] = %.2f\n", i, arr[i]);
	}

	//printf("x = %.2f\ny = %.2f\nz = %.2f\n", x, y, z);

	return 0;
}