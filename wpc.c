#include <math.h>
#include "stdio.h"
int validation_mode = 1; // 启用验证

int **swap_row(int **mat, int row1, int row2, int col) // 互换row1,row2行
{
	int i, tmp;
	for (i = 0; i < col; i++)
	{
		tmp = mat[row1][i];
		mat[row1][i] = mat[row2][i];
		mat[row2][i] = tmp;
	}
	return mat;
}

void up_triangle(int **mat, int row, int col)
{
	int i, j, k;
	if (row <= 1) return;
	k = 0; // 1的数量
	for (j = 0; j < col; j++)
	{
		for (i = 0; i < row; i++)
		{
			if(mat[i][j]) mat = swap_row(mat, i, k++, col);// 互换i行和k行
		}
		if (k) break;
	}
	// 化简后k-1行
	for (i = 1; i < k; i++)
		for (j = 0; j < col; j++)
			mat[i][j] = abs((mat[i][j] - mat[0][j]) % 2);
    
	up_triangle(&mat[1], row - 1, col);
}

int down_triangle(int **mat, int row, int col)
{
	int i, j, k;
	for (k = 0; k < row - 1; k++)
		for (i = k + 1; i < row; i++)
			for (j = 0; j < col; j++)
				if (mat[i][j])
				{
					if (mat[k][j])
						for (; j < col; j++)
							mat[k][j] = abs((mat[k][j] - mat[i][j]) % 2);
					break;
				}
	return 1;
}

int valid_up_triangle(int **mat, int row, int col) // mat=(A,b), i=rank(A), j=rank(b)
{
	int i, j;
	for (i = 0, j = 0; j < col - 1 && i < row; j++)
	{
		if (mat[i][j]) i++;
	}
	for (j = row; j > 0; j--) if (mat[j - 1][col - 1]) break;
	if ((!validation_mode && i < j) || (validation_mode && i < row)) return 0;
	else return 1;
}

int *mat_solve(int **A, int *b, int row, int col) //求解Ax=b
{
	int i, j,
	    **Ab = (int**)calloc(row, sizeof(int*)), // 增广矩阵
	    *solution = (int*)calloc(col, sizeof(int)); // Ax=b的解

	// 构建增广矩阵Ab
	for (i = 0; i < row; i++)
	{
		Ab[i] = (int*)calloc(col + 1, sizeof(int));
		for (j = 0; j < col; j++)
		{
			Ab[i][j] = A[i][j];
		}
		Ab[i][col] = b[i];
	}
	// 求解Ax=b
	up_triangle(Ab, row, col + 1);
	if(!valid_up_triangle(Ab, row, col + 1)) // 无解
	{
		solution[0] = -1;
		return solution;
	}
	down_triangle(Ab, row, col + 1);

	for (i = 0; i < row; i++)
	{
		for (j = i; j < col; j++)
			if (Ab[i][j])
			{
				solution[j] = Ab[i][col];
				break;
			}
	}
	return solution;
}

int WPC_encode_block(int **H, int *cover, int *stego, int *wet_indicate, int blocklen, int *msg, int *msglen)
{
	int i, j, k, wet_num = 0, msg_len = *msglen;
	for (i = 0; i < blocklen; i++)
		if (wet_indicate[i]) wet_num++;

	int **H_wet = (int**)calloc(msg_len, sizeof(int*)), // 分段湿校验矩阵
		**H_dry = (int**)calloc(msg_len, sizeof(int*)), // 分段干校验矩阵
		*msg_wet = (int*)calloc(msg_len, sizeof(int)),
		*msg_dry = (int*)calloc(msg_len, sizeof(int)),
		*cover_wet = (int*)calloc(wet_num, sizeof(int)),
		*stego_dry;

	for (i = 0; i < msg_len; i++)
	{
		H_wet[i] = (int*)calloc(wet_num, sizeof(int));
		H_dry[i] = (int*)calloc(blocklen - wet_num, sizeof(int));
		for (j = 0, k = 0; j < blocklen; j++)
		{
			if (wet_indicate[j]) H_wet[i][k++] = H[i][j];
			else H_dry[i][j - k] = H[i][j];
		}
	}
	for (j = 0, k = 0; j < blocklen; j++) if (wet_indicate[j]) cover_wet[k++] = cover[j];

	for (i = 0; i < msg_len; i++) // msg_wet = H_wet * cover_wet
	{
		for (j = 0; j < wet_num; j++)
			msg_wet[i] += H_wet[i][j] * cover_wet[j];
		msg_wet[i] %= 2;
		msg_dry[i] = msg_wet[i] == msg[i] ? 0 : 1;
	}
	
	stego_dry = mat_solve(H_dry, msg_dry, msg_len, blocklen - wet_num);
	if (stego_dry[0] == -1) // 无解
	{
		if(validation_mode) *msglen = 0;
		for (i = 0; i < blocklen; i++) stego[i] = cover[i];
	}
	else
		for (i = 0, k = 0; i < blocklen; i++) stego[i] = wet_indicate[i] ? cover[i] : stego_dry[k++];
	
	for (i = 0; i < msg_len; i++) free(H_wet[i]);
	for (i = 0; i < msg_len; i++) free(H_dry[i]);
	free(H_wet);
	free(H_dry);
	free(msg_wet);
	free(msg_dry);
	free(cover_wet);
	free(stego_dry);
	return 0;
}

int* WPC_encode(int **H, int blocklen, int *cover, int *wet_indicate, int coverlen, int *msg, int msglen, int *msgembedlen)
{
	int i, j, blocknum = floor(coverlen / blocklen), drynum = 0, msgpos = 0, msgblklen,
		*stego = (int*)calloc(coverlen, sizeof(int)),
		*stego_blk = (int*)calloc(blocklen, sizeof(int));

	for (i = 0; i < coverlen; i++) drynum += 1 - wet_indicate[i];

	for (i = 0; i < blocknum; i++)
	{
		for (msgblklen = j = 0; j < blocklen; j++) msgblklen += 1 - wet_indicate[i*blocklen+j];
		msgblklen = round(msgblklen * 1.0 / drynum * msglen); // 载荷分配
		WPC_encode_block(H, &cover[i*blocklen], stego_blk, &wet_indicate[i*blocklen], blocklen, &msg[msgpos], &msgblklen);
		memcpy(&stego[i*blocklen], stego_blk, blocklen * sizeof(int));
		msgpos += msgblklen;

		if (msgpos > msglen)
		{
			memcpy(&stego[i*blocklen + blocklen], &cover[i*blocklen + blocklen], (coverlen - i * blocklen - blocklen) * sizeof(int));
			break;
		}
	}
	*msgembedlen = msgpos;
	printf("embed %d of %d bits.\n", msgpos, msglen);
	free(stego_blk);
	return stego;
}

int WPC_decode_block(int **H, int *stego, int *msg, int msglen, int blocklen)
{
	int i, j;
	for (i = 0; i < msglen; i++) // msg = H * stego
	{
		msg[i] = 0;
		for (j = 0; j < blocklen; j++)
			msg[i] += H[i][j] * stego[j];
		msg[i] %= 2;
	}
}

int* WPC_decode(int **H, int blocklen, int *stego, int *wet_indicate, int coverlen, int msglen)
{
	int i, j, blocknum = floor(coverlen / blocklen), drynum = 0, msgpos = 0, drynumblk, msgblklen,
		*msg = (int*)calloc(msglen, sizeof(int)),
		*msg_blk = (int*)calloc(blocklen, sizeof(int)), **H_valid;

	if (validation_mode)
	{
		H_valid = (int**)calloc(blocklen, sizeof(int*));
		for (i = 0; i < blocklen; i++)
			H_valid[i] = (int*)calloc(blocklen+1, sizeof(int));
	}

	for (i = 0; i < coverlen; i++) drynum += 1 - wet_indicate[i];

	for (i = 0; i < blocknum; i++)
	{
		for (drynumblk = j = 0; j < blocklen; j++) drynumblk += 1 - wet_indicate[i*blocklen + j]; // dry point num
		msgblklen = round(drynumblk * 1.0 / drynum * msglen); // 载荷分配
		if(validation_mode)
		{
			for (j = 0; j < msgblklen; j++)
				for (int k1 = 0, k2 = 0; k1 < blocklen; k1++)
					if (!wet_indicate[i*blocklen + k1]) H_valid[j][k2++] = H[j][k1];
			up_triangle(H_valid, msgblklen, drynumblk + 1);
			if (!valid_up_triangle(H_valid, msgblklen, drynumblk + 1)) msgblklen = 0; // 无解
		}
		if (msgpos + msgblklen > msglen) break;
		WPC_decode_block(H, &stego[i*blocklen], msg_blk, msgblklen, blocklen);
		memcpy(&msg[msgpos], msg_blk, msgblklen * sizeof(int));
		msgpos += msgblklen;
	}

	free(msg_blk);
	return msg;
}

int main()
{
	int i, j, blocklen = 4, coverlength = 40000, msglen = 8000, msgembedlen, dec_err,
		**H = (int**)calloc(blocklen, sizeof(int*)), // 校验矩阵
	    *cover = (int*)calloc(coverlength, sizeof(int)),
		*wet_indicate = (int*)calloc(coverlength, sizeof(int)),
		*msg = (int*)calloc(msglen, sizeof(int)),
		*stego, *msg_dec;
	
	srand(time(0));
	for (i = 0; i < blocklen; i++)
	{
		H[i] = (int*)calloc(blocklen, sizeof(int));
		for (j = 0; j < blocklen; j++)
			H[i][j] = rand() % 2;
	}
	for (i = 0; i < coverlength; i++)
	{
		if(i < msglen) msg[i] = rand() % 2;
		cover[i] = rand() % 2;
		wet_indicate[i] = rand() % 2;
	}

	stego = WPC_encode(H, blocklen, cover, wet_indicate, coverlength, msg, msglen, &msgembedlen);
	msg_dec = WPC_decode(H, blocklen, stego, wet_indicate, coverlength, msglen);
	for (dec_err = 0, i = 0; i < msgembedlen; i++) dec_err += msg_dec[i] == msg[i] ? 0 : 1;
	printf("decode error num: %d\n", dec_err);

	for (i = 0; i < blocklen; i++) free(H[i]);
	free(H);
	free(cover);
	free(stego);
	free(wet_indicate);
	free(msg);
	free(msg_dec);
	return 0;
}