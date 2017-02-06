/*

 * 2DPotts_parallel.cpp
 *
 *  Created on: 2013年10月25日
 *      Author: BS
 */
//该程序目的是将2DPoots并行化
//本程序做标准Potts模型
/*
 2013/5/5
 16:59
 re-write function makeup() by Woleff

 2013/11/12
 7:17
 parallel program is right
 But need to add function that can insert temperature
 we can put it in the "wan" function.
 make an array save the insert temperature.

 2013/12/22
 overlap函数仍然存在问题

 2013/12/29
 研究overlap函数
 21:40
 overlap函数应该没有问题了
 下面需要写F_P函数

 2014/1/11
 0:00
 态密度g计算存在问题

 2014/1/12
 20:54
 出现每个温度不独立的情况，导致无法并行运算

 2014/1/18
 22;16
 内存与m有关，这种情况不应该发生
 */

#include<iostream>
#include <fstream>
#include <time.h>
#include<fstream>
#include<iomanip>
#include<math.h>
#include<stdlib.h>
#include <sstream>
#include <omp.h>

using namespace std;

const int num_E = 10000;

struct point //每个个点上的点
{
	int num;
	int x; //储存数组坐标 从0开始
	int y;
	bool cc; //代表是否被makeup储存 储存为1 为储存为0
};

int *seed = new int[omp_get_num_procs()];
void begin_seed(int * seed_begin){
	for (int i = 0; i < omp_get_num_procs(); i++){
		seed_begin[i] = i + 1;
	}
}

const int N = 10; //为二维数组的长度和宽度
const double J = 1;
const int nn = 1; const int q = 2; //potts模型取值范围

struct EEx {
	double EE;
	double gE;
};

/*
int main(){
begin_seed(seed);
void F_P(int Nequi, int m, int N0, double T, double B, double P, EEx *EP_pre_T,
int &num_EP_prese_T);
int N0_sure(double T);

int Nequi =1000;
int m = 30000;
int N0 = 0;
double T_i[13] = {0.18,0.24,0.42,0.48,0.54,0.6,0.72,0.9,1.02,1.08,1.2,1.32,1.5};

#pragma omp parallel for
for(int i=0;i<13;i++){
EEx * EP = new EEx [1000];
int num_EP = 0;

double T = T_i[i];
double B = 1.0/T;
double P = 1.0 - exp(-B * J);
N0 = N0_sure(T);

F_P(Nequi,m,N0,T,B,P,EP,num_EP);
stringstream sT;
sT << T;

ofstream fout(("EP_T_" + sT.str() + "_N=32.txt").c_str());
for(int j=0;j<num_EP;j++){
fout<<EP[j].EE/double(N*N)<<"\t"<<setprecision(16)<<EP[j].gE/double(m)<<endl;
}

fout.close();
delete [] EP;
}

return 0;
}
*/
//F-Sigma2
/*
int main() //打算做q=2的potts模型的E值
{
begin_seed(seed);
void F_sigma(int Nequi, double T_start, double T_end, int xunhuan, int num_T,
int Sigma_i_T,int num_m);

omp_set_nested(1);

int Nequi = 1000;

int num_T = 75;
F_sigma(Nequi,0,6,32,num_T,num_T-1,10);

return 0;
}
*/





double abs_local(double x) {
	return x >= 0 ? x : -1.0 * x;
}

int N0_sure(double T) {

	//return 3;
	if (T <= 0.2)
		return 1;
	else if (T > 6) {
		return 26;
	}
	else {
		if (4.0 * T - int(4.0 * T) >= 0.5)
			return int(4.0 * T) + 1;
		else
			return int(4.0 * T);
	}
}

int myround(double x) {
	//double add=0.5;
	//int min,max;
	int sa;
	double si;	    //sa用于保存x的整数部分,si用于保存加0.5后的临时值。整数部分的四舍五入
	if (x == 0.0)
		return 0;
	else if (x > 0.0) {
		sa = (int)x;
		si = x + 0.5;
		if (sa == floor(si))	    //如果if语句成立说明x的小数比0.5小，应当舍去
			return sa;
		else
			return sa + 1;
	}
	else  //负数部分的四舍五入
	{
		sa = (int)x;
		si = x - 0.5;
		if (sa == ceil(si))
			return sa;
		else
			return sa - 1;
	}
}

int rand_16807(int & seed_local_rand)
{
	const int a = 16807;         // 16807 法
	const int b = 0;
	const int m = 2147483647;    // MAX_INT
	const int q = 127773;        // q = m / a;
	const int r = 2836;          // r = m % a;
	int _z = a * (seed_local_rand % q) - r * (int)(seed_local_rand / q) + b; // 计算 mod
	if (_z < 0)
		_z += m;             // 将结果调整到 0 ~ m
	return seed_local_rand = _z;
}

double randw() {


	return (double(rand_16807(seed[omp_get_thread_num()])) / 2147483647);


	//return (double(rand()) / RAND_MAX);
}

//产生随机数0,1
int randm() {
	int r;

	r = randw() >= 0.5 ? 1 : 0;

	return r;
}

//产生随机数1-q           //q在最开始声明
int randN(int N_local) {
	int r;

	r = myround(double(N_local)*randw() - 0.5);

	if (r == N_local)
		return N_local - 1;
	if (r == -1)
		return 0;

	return r;
}

//2->10
int zh2A(point s[], int length_s) {
	int value = 0;
	int w = 1;
	for (int i = 0; i < length_s; i++) {
		value +=
			(s[length_s - i - 1].num == 0) ?
			0 : s[length_s - i - 1].num * w;
		w *= 2;
	}

	return value;
}

void average(double s[], int length, double &value_average, double& sted) {
	value_average = 0;
	sted = 0;
	for (int i = 0; i < length; i++) {
		value_average += s[i] / double(length);
		sted += s[i] / double(length) * s[i];
	}
	sted = sqrt(sted - value_average * value_average);

}

//横向输出EP
void printEP(EEx **EP, int *num_EP, int num_T, int num_Ptxt, int m){
	//找到num最大值
	int num_max = 0;
	for (int i = 0; i < num_T; i++){
		if (num_max < num_EP[i])
			num_max = num_EP[i];
	}

	stringstream stxt;
	stxt << num_Ptxt;
	//输出
	ofstream fout_p(("p" + stxt.str() + ".txt").c_str());
	for (int i = 0; i < num_max; i++){
		for (int j = 0; j<num_T; j++){
			if (j != num_T - 1){
				if (i>num_EP[j]){
					fout_p << 0 << "\t" << 0 << "\t";
				}
				else if (i == num_EP[j]){
					fout_p << EP[j][i].EE << "\t" << 0 << "\t";
				}
				else{
					if (EP[j][i].gE>0)
						fout_p << EP[j][i].EE << "\t" << EP[j][i].gE / double(m) << "\t";
					else
						fout_p << 0 << "\t" << 0 << "\t";
				}
			}
			else{
				if (i > num_EP[j]){
					fout_p << 0 << "\t" << 0;
				}
				else if (i == num_EP[j]){
					fout_p << EP[j][i].EE << "\t" << 0;
				}
				else{
					if (EP[j][i].gE > 0)
						fout_p << EP[j][i].EE << "\t" << EP[j][i].gE / double(m);
					else
						fout_p << 0 << "\t" << 0;
				}
			}

		}
		fout_p << endl;
	}
}

void printEP2(EEx **EP, int *num_EP, int num_T_total){
	ofstream fout_p3("p3.txt");
	for (int i = 0; i < num_T_total; i++){
		for (int j = 0; j < num_EP[i]; j++){
			fout_p3 << i << "\t" << EP[i][j].EE << "\t" << EP[i][j].gE << endl;
		}
	}
}

//赋初值
void begin(point **s) {
	//cout<<"begin"<<endl;

	for (int bi = 0; bi < N; bi++) {
		for (int bj = 0; bj < N; bj++) {
			s[bi][bj].num = randN(q);
			s[bi][bj].x = bi;
			s[bi][bj].y = bj;
			s[bi][bj].cc = false;
		}
	}

}

void makeup(point **team, double P) {


	point *upstore[N * N];
	int len_save_store_ch = 0;  //save the point have changed
	int len_save_store_tt = 0;  //save the poind choised but not changed
	int dx, dxn, dxa, dy, dya, dyn;

	for (int mi = 0; mi < N; mi++) {
		for (int mj = 0; mj < N; mj++)
			team[mi][mj].cc = false;
	}

	upstore[0] = &team[randN(N)][randN(N)];
	(*upstore[0]).cc = true;
	len_save_store_tt++;
	int qx = randN(q);
	do {


		dx = (*upstore[len_save_store_ch]).x;
		dy = (*upstore[len_save_store_ch]).y;
		dxa = (dx == N - 1) ? 0 : (dx + 1);
		dya = (dy == N - 1) ? 0 : (dy + 1);
		dxn = (dx == 0) ? N - 1 : (dx - 1);
		dyn = (dy == 0) ? N - 1 : (dy - 1);


		if ((team[dx][dya].cc == false) && (randw() < P)
			&& ((*upstore[len_save_store_ch]).num == team[dx][dya].num)) //向右检查
		{

			upstore[len_save_store_tt] = &team[dx][dya];
			team[dx][dya].cc = true;
			len_save_store_tt++;

		}


		if ((team[dxa][dy].cc == false) && (randw() < P)
			&& ((*upstore[len_save_store_ch]).num == team[dxa][dy].num)) //向下检查
		{

			upstore[len_save_store_tt] = &team[dxa][dy];
			team[dxa][dy].cc = true;
			len_save_store_tt++;

		}


		if ((team[dxn][dy].cc == false) && (randw() < P)
			&& ((*upstore[len_save_store_ch]).num == team[dxn][dy].num)) //向上检查
		{

			upstore[len_save_store_tt] = &team[dxn][dy];
			team[dxn][dy].cc = true;
			len_save_store_tt++;

		}


		if ((team[dx][dyn].cc == false) && (randw() < P)
			&& ((*upstore[len_save_store_ch]).num == team[dx][dyn].num)) //向左检查
		{

			upstore[len_save_store_tt] = &team[dx][dyn];
			team[dx][dyn].cc = true;
			len_save_store_tt++;

		}

		(*upstore[len_save_store_ch]).num = qx;
		len_save_store_ch++;   //
		if (len_save_store_ch == len_save_store_tt) {
			break;
		}
	} while (1);

}

void jump(int n, point **team, double P) {

	for (int ji = 0; ji < n; ji++) {
		makeup(team, P);
	}

}

//周期循环
double H(point **s) {
	double HH1 = 0;

	for (int hi = 0; hi < N; hi++) //先算横的
	{
		for (int hj = 0; hj < N - 1; hj++) {
			HH1 += double(s[hi][hj].num) == double(s[hi][hj + 1].num) ? 1.0 : 0;
		}
		HH1 += double(s[hi][N - 1].num) == double(s[hi][0].num) ? 1.0 : 0;
	}

	double HH2 = 0;

	for (int hj = 0; hj < N; hj++) //纵的
	{
		for (int hi = 0; hi < N - 1; hi++) {
			HH2 += double(s[hi][hj].num) == double(s[hi + 1][hj].num) ? 1.0 : 0;
		}
		HH2 += double(s[N - 1][hj].num) == double(s[0][hj].num) ? 1.0 : 0;
	}

	return -1.0 * (J * (HH1 + HH2));
}

double get_T(double T_start, double T_end, int i_T, int num_T_total) {
	double value_T = 0;

	value_T = (T_end - T_start) / double(num_T_total) * (i_T + 1.0) + T_start;

	return value_T;
}

double get_T_high_temperature(double T_start, double T_end, int i_T, int num_T_total) {
	double value_T = 0;

	if (i_T != num_T_total - 1){
		value_T = (T_end - T_start) / double(num_T_total - 1) * (i_T + 1.0) + T_start;
	}
	else{
		value_T = 1.E10;
	}

	return value_T;
}


double comper(EEx *EP1, int &num_EP1, EEx *EP2, int &num_EP2, int m){
	//主要想法：对于一个类抛物线分布，需要找到顶点。

	//1.找顶点:利用加权平均的方法，sigma:P*EE,求平均，认为就是顶点{---------------------
	double topE1 = 0;
	double topE2 = 0;


	topE1 = EP1[num_EP1].EE;

	topE2 = EP2[num_EP2].EE;
	//end1}=======================================================================

	//2.寻找交点.利用差值交换符号{---------------------------------------------------

	//利用EP_left 或right 代替原有{
	EEx * EP_left;
	int *num_left;
	double topE_left;
	EEx * EP_right;
	int *num_right;
	double topE_right;
	if (topE1 <= topE2){
		EP_left = EP1;
		num_left = &num_EP1;
		topE_left = topE1;
		EP_right = EP2;
		num_right = &num_EP2;
		topE_right = topE2;
	}
	else{
		EP_left = EP2;
		num_left = &num_EP2;
		topE_left = topE2;
		EP_right = EP1;
		num_right = &num_EP1;
		topE_right = topE1;
	}

	for (int i = 0; i < *num_left; i++){
		if (topE_left <= EP_left[i].EE){
			if (i != 0)
				topE_left = EP_left[i - 1].EE;

			break;
		}
	}

	for (int i = 0; i < *num_right; i++){
		if (topE_right <= EP_right[i].EE){
			topE_right = EP_right[i].EE;
			break;
		}
	}
	//利用EP_left 或right 代替原有}==================================

	//左支的左半面一定保留，但对于右支则选择交点之后保留，若始终没有交点，则右支可以删除

	int i_mid_left = -1;//中点在左边的数
	int i_mid_right = -1;//中点在右边的数

	if (EP_left[*num_left - 1].EE >= EP_right[0].EE){
		int i = 0;
		for (; i < *num_left; i++){

			double E_local = EP_left[i].EE;
			if (E_local >= topE_left){
				int i_E_local_right = 0;
				for (; i_E_local_right < *num_right; i_E_local_right++){
					if (E_local == EP_right[i_E_local_right].EE){//&& EP_right[i_E_local_right].gE > 0){
						break;
					}
				}

				if (i_E_local_right == *num_right)
					continue;

				if (fabs(EP_left[i].gE) <= fabs(EP_right[i_E_local_right].gE)){
					if (i == 0 || i_E_local_right == 0){
						i_mid_left = i;
						i_mid_right = i_E_local_right;
					}
					else{

						//选择交点，是该点还是变化之前的点
						if (EP_left[i - 1].EE == EP_right[i_E_local_right - 1].EE){
							//&& EP_left[i-1].gE > 0 && EP_right[i_E_local_right-1].gE > 0){
							double P1_left = fabs(EP_left[i].gE);
							double P1_right = fabs(EP_right[i_E_local_right].gE);
							double P2_left = fabs(EP_left[i - 1].gE);
							double P2_right = fabs(EP_right[i_E_local_right - 1].gE);
							double deta1_local = fabs(P1_left - P1_right);
							double deta2_local = fabs(P2_left - P2_right);

							bool b_1 = 1;
							//deta + 1 防止deta为0
							double pd_1 = P1_left * P1_right / (deta1_local + 0.1);
							double pd_2 = P2_left * P2_right / (deta2_local + 0.1);
							if (pd_1 >= pd_2) {
								b_1 = true;
							}
							else {
								b_1 = false;
							}

							if (b_1){
								i_mid_left = i;
								i_mid_right = i_E_local_right;
							}
							else{
								i_mid_left = i - 1;
								i_mid_right = i_E_local_right - 1;
							}
						}

						else{
							i_mid_left = i;
							i_mid_right = i_E_local_right;
						}

					}
				}
				/*
				else if(i == (*num_left-1) || EP_left[i+1].gE <0){
				i_mid_left = i;
				i_mid_right = i_E_local_right;
				}
				*/
			}
			if (i_mid_left != -1)
				break;
		}
		if (i == *num_left){//说明右支始终在左支下面
			if (EP_left[*num_left - 1].EE >= EP_right[*num_right - 1].EE){
				//完全删除右支
				i_mid_right = 0;
				i_mid_left = *num_left - 1;
			}
			else{
				i_mid_left = *num_left;
				do
				{
					i_mid_left--;
					for (int j = 0; j < *num_right; j++){
						if (EP_left[i_mid_left].EE == EP_right[j].EE){
							i_mid_right = j;
						}
					}
				} while (i_mid_right == -1 && i_mid_left >= 0);
			}

		}
	}
	//end2}==============================================


	//3.将交点以外的点的概率取负，这样以此判断{-----------------------------

	if (i_mid_right != -1){
		for (int i = i_mid_left + 1; i < (*num_left); i++){
			if (EP_left[i].gE > 0){

				for (int j = i_mid_right; j < *num_right; j++){//只有当另一条含有该条的点时才应该删除

					//if(EP_right[j].EE == EP_left[i].EE)
						{

							//if( EP_right[j].gE>0)
							{
								EP_left[i].gE *= -1;
								break;
							}
						}//else
					if (EP_right[j].EE > EP_left[i].EE){
						break;
					}

				}
				/*
				if(b_de)
				EP_left[i].gE *= -1;
				*/
			}
		}
		for (int i = 0; i < i_mid_right; i++){
			if (EP_right[i].gE > 0){

				for (int j = 0; j < i_mid_left; j++){
					//if(EP_left[j].EE == EP_right[i].EE )
						{

							//if (EP_left[j].gE > 0)
							{

								EP_right[i].gE *= -1;
								break;
							}
						}//else
					if (EP_left[j].EE > EP_right[i].EE){
						break;
					}
					/*
					if(b_de )
					EP_right[i].gE *= -1;
					*/
				}
			}
		}
	}



	/*
		cout<<"........................."<<endl;
		cout<<"left_A\t"<<EP_left[i_mid_left].EE<<endl;
		for(int i=0;i<*num_left;i++){
		cout<<EP_left[i].EE<<"\t"<<EP_left[i].gE<<endl;
		}
		cout<<"\nright_A\t"<<EP_right[i_mid_right].EE<<endl;
		for(int i=0;i<*num_right;i++){
		cout<<EP_right[i].EE<<"\t"<<EP_right[i].gE<<endl;
		}
		cout<<"---------------------"<<endl;
		*/
	return 0;
}

bool CheckinOtherT(EEx **EP, int num_T, int *num_EP, int i_T){
	/*
	 * 该函数用于检查i_T温度的EP曲线两端点是否和其他曲线存在交点
	 */

	//首先找到曲线头和尾
	int S = -1;
	int e = -1;
	int num_EPoverzero = 0;
	for (int i = 0; i < num_EP[i_T]; i++){
		if (EP[i_T][i].gE > 0){
			num_EPoverzero++;
			if (i == 0){
				S = 0;
			}
			else if (EP[i_T][i - 1].gE < 0){
				S = i;
			}
			if (i == num_EP[i_T] - 1){
				e = num_EP[i_T] - 1;
			}
			else if (EP[i_T][i + 1].gE < 0){
				e = i;
			}
		}
		if (S != -1 && e != -1){
			break;
		}
	}

	if (num_EPoverzero == 0)
		return true;



	//寻找头和尾是否在别的温度中{------------------------------------------------
	bool b_findstart = false;
	bool b_findend = false;
	//两个特殊情况，位于整个曲线的头和尾是不可能有overlap的

	if (EP[i_T][S].EE == -2 * N*N){
		b_findstart = true;
	}
	double Emax = -2.0*N*N;
	for (int i_check = 0; (i_check < num_T); i_check++){
		if (i_check == i_T)
			continue;

		for (int j_check = 0; j_check < num_EP[i_check]; j_check++){
			if (EP[i_check][j_check].gE > 0){
				if (EP[i_check][j_check].EE > Emax){
					Emax = EP[i_check][j_check].EE;
				}
				if (EP[i_check][j_check].EE == EP[i_T][S].EE)
					b_findstart = true;
				if (EP[i_check][j_check].EE == EP[i_T][e].EE)
					b_findend = true;
			}
			if (b_findend && b_findstart)
				break;
		}
		if (b_findend && b_findstart)
			break;
	}

	if (Emax <= EP[i_T][e].EE){
		b_findend = true;
	}



	//寻找overlapend}=====================================================

	//如果没有找到头和尾,需要增添{-----------------------------------------------
	/*
		if(!b_findstart){
		for(int i=0;i<S;i++){
		}
		if(S == 0){
		cout<<"erro:S=0,there is no overlap before"<<endl;
		exit(0);
		}

		}
		if(!b_findend){
		if(e == (num_EP[i_T]-1)){
		cout<<"erro:e=num,there is no overlap after"<<endl;
		exit(0);
		}

		}
		*/
	if (b_findend && b_findstart){
		return true;
	}
	else{
		EP[i_T][(num_EP[i_T] + 1)].EE = S;
		EP[i_T][(num_EP[i_T] + 2)].EE = e;
		if (!b_findend)
			EP[i_T][(num_EP[i_T] + 3)].EE = 1.0;
		else
			EP[i_T][(num_EP[i_T] + 3)].EE = 0.0;
		if (!b_findstart)
			EP[i_T][(num_EP[i_T] + 4)].EE = 1.0;
		else
			EP[i_T][(num_EP[i_T] + 4)].EE = 0.0;
		return false;
	}


}

//单独一个函数，作用是找到overlap，在overlap一下的能量将淘汰掉

void overlap(EEx **EP, int *num_EP, int num_T_total, int m){


	//倆倆比較
	for (int i = 0; i < (num_T_total - 1); i++){
		for (int j = i + 1; (j < num_T_total); j++){
			/*
			if((i == 61 || j ==61 )&& i>= 16){
			printEP(EP,num_EP,num_T_total,1);
			}
			*/
			comper(EP[i], num_EP[i], EP[j], num_EP[j], m);


			//查看是否只有一个大于零，将其删除
			int num_EP_local = 0;
			int location_EP = 0;
			for (int i_cp = 0; i_cp < num_EP[i]; i_cp++){
				if (EP[i][i_cp].gE > 0){
					num_EP_local++;
					location_EP = i_cp;
				}
			}
			if (num_EP_local == 1)
				EP[i][location_EP].gE *= -1;

			num_EP_local = 0;
			for (int i_cp = 0; i_cp < num_EP[j]; i_cp++){
				if (EP[j][i_cp].gE > 0){
					num_EP_local++;
					location_EP = i_cp;
				}
			}
			if (num_EP_local == 1)
				EP[j][location_EP].gE *= -1;
			/*
			if(i == 61 || j ==61 ){
			printEP(EP,num_EP,num_T_total,2);
			}
			*/

		}
		//printEP(EP,num_EP,num_T_total,3,m);
	}
	/*
		ofstream fout_p3("p3.txt");
		for(int i=0;i<num_T_total;i++){
		for(int j=0;j<num_EP[i];j++){
		fout_p3<<i<<"\t"<<EP[i][j].EE<<"\t"<<EP[i][j].gE<<endl;
		}
		}
		*/
	//	printEP(EP,num_EP,num_T_total,2);
	//检查是否断开{----------------------------------------------------------------
	//由于覆盖方式采取两两比较的方法，所以可以对每一段温度进行检查，检查出现P分布曲线两端点是否在别的曲线中
	int arry_check[300];
	int leng_arry_check = 0;
	for (int i = 0; i < 300; i++){
		arry_check[i] = -1;
	}
	for (int i_T = 0; i_T < num_T_total; i_T++){
		if (!(CheckinOtherT(EP, num_T_total, num_EP, i_T))){
			arry_check[leng_arry_check] = i_T;
			leng_arry_check++;
		}
	}
	//end}======================================================================
	//将断开的重新不上{--------------------------------------------------------------

	do{
		int i_T1 = 0;
		int i_T2 = 0;
		int positionT1inArry = 0;
		int positionT2inArry = 0;
		int i_T1_e = 0;
		int i_T2_s = 0;
		double Ee = 0;
		for (int i = 0; i<leng_arry_check; i++){
			int i_T = arry_check[i];
			if (EP[i_T][num_EP[i_T] + 3].EE == 1.0){
				if (Ee > EP[i_T][myround(EP[i_T][num_EP[i_T] + 2].EE)].EE){//change int
					Ee = EP[i_T][myround(EP[i_T][num_EP[i_T] + 2].EE)].EE;
					i_T1 = i_T;
					i_T1_e = myround(EP[i_T][num_EP[i_T] + 2].EE);
					positionT1inArry = i;
				}
			}
		}
		Ee = 0;
		for (int i = 0; i<leng_arry_check; i++){
			int i_T = arry_check[i];
			if (EP[i_T][num_EP[i_T] + 4].EE == 1.0){
				if (Ee > EP[i_T][myround(EP[i_T][num_EP[i_T] + 1].EE)].EE && Ee > EP[i_T1][i_T1_e].EE){
					Ee = EP[i_T][myround(EP[i_T][num_EP[i_T] + 1].EE)].EE;
					i_T2 = i_T;
					i_T2_s = myround(EP[i_T][num_EP[i_T] + 1].EE);
					positionT2inArry = i;
				}
			}
		}



		//尝试一下直接连接
		//	fflush(stdout);
		//	cout<<i_T1<<"\t"<<i_T2<<endl;
		//	fflush(stdout);
		for (int i = i_T1_e + 1; i < num_EP[i_T1]; i++){
			if (EP[i_T1][i].gE < 0){
				EP[i_T1][i].gE *= -1;
			}
		}
		for (int i = 0; i < i_T2_s; i++){
			if (EP[i_T2][i].gE < 0){
				EP[i_T2][i].gE *= -1;
			}
		}
		comper(EP[i_T1], num_EP[i_T1], EP[i_T2], num_EP[i_T2], m);



		//删除掉两个点
		if (EP[i_T1][num_EP[i_T1] + 4].EE != 1){
			leng_arry_check--;
			for (int i = positionT1inArry; i < (leng_arry_check); i++){
				arry_check[i] = arry_check[i + 1];
			}
		}
		else{
			EP[i_T1][num_EP[i_T1] + 3].EE = 0;
		}

		if (positionT1inArry < positionT2inArry)
			positionT2inArry--;
		if (EP[i_T2][num_EP[i_T2] + 3].EE != 1){
			leng_arry_check--;

			for (int i = positionT2inArry; i < (leng_arry_check); i++){
				arry_check[i] = arry_check[i + 1];
			}
		}
		else{
			EP[i_T2][num_EP[i_T2] + 4].EE = 0;
		}

		//	printEP(EP,num_EP,num_T_total,2);
	} while (leng_arry_check > 0);
}


//nn为数组位数 EP为储存能量E和不同温度下概率P num_T从0开始
void F_P(int Nequi, int m, int N0, double T, double B, double P, EEx *EP_pre_T,
	int &num_EP_prese_T) {


	void BubbleSort(EEx * pData, int Count);

	point **team;
	team = new point *[N];
	for (int i = 0; i < N; i++) {
		team[i] = new point[N];
	}
	int leng = 0;

	EEx *EP_local = new EEx[num_E];

	//开始
	double HH = 0;

	begin(team);

	jump(Nequi, team, P);

	for (int ii = 0; ii < m; ii++) {


		jump(N0, team, P);

		//计算
		HH = H(team);
		int i = 0;
		for (; i < leng; i++) {
			if (EP_local[i].EE == HH)
				break;
		}

		if (i == leng) {
			EP_local[leng].EE = HH;
			EP_local[leng].gE = 1;
			leng++;
			if (leng >= num_E){
				printf("erro: num_E not enough.");
				exit(0);
			}
		}
		else {
			EP_local[i].gE++;
		}

	}
	//排序{-------------------------------------------

	BubbleSort(EP_local, leng);

	//排序end}========================================


	//找到最大能量{------------------------------------------------------

	double topE = 0;
	double topEgE = 0;

	for (int i_top = 0; i_top < leng; i_top++){
		if (topEgE < EP_local[i_top].gE){
			topE = EP_local[i_top].EE;
			topEgE = EP_local[i_top].gE;
		}
	}
	EP_local[leng].EE = topE;

	//找到最大能量}=========================================================

	//将局部变量复制到传递进来的变量{--------------------------------------------
	num_EP_prese_T = 0;

	for (int i = 0; i < leng; i++) {
		EP_pre_T[i].EE = EP_local[i].EE;
		EP_pre_T[i].gE = EP_local[i].gE;
	}
	num_EP_prese_T = leng;
	EP_pre_T[num_EP_prese_T].EE = EP_local[leng].EE;

	if (num_EP_prese_T == 0){
		printf("~~~~~~~~~~erro~~~~~~~~~~~\n");
		printf("~~~~~~~~~~~~~~~~~~~~~~~~~\n");
		fflush(stdout);
		exit(0);
	}
	//end}======================================================================
	for (int i = 0; i < N; i++)
		delete[] team[i];
	delete[] team;

	delete[] EP_local;
}

//将原all函数分开 num_T_total是总数需要加一
void P_F(EEx **EP, int *num_EP, int num_T_total, int Nequi, int m,
	double T_start, double T_end, double **FNNf, int xunhuan) {

	const int leng_ENE = 2 * N * N + 10;
	EEx EE[leng_ENE];
	int leng_EE = 0;
	//初始化

	for (int i = 0; i < leng_ENE; i++) {
		EE[i].EE = 1;
		EE[i].gE = 0;
	}
	leng_EE = 0;

	//step 1 overlap 将曲线连接起来
	//step 2 计算g
	//step 3 计算f


	//overlap{-----------------------------------------------------------------
	//	printEP(EP,num_EP,num_T_total,1);

	overlap(EP, num_EP, num_T_total, m);



	//	printEP(EP,num_EP,num_T_total,2);
	//overlap}====================================================================



	//计算g{========================================================================

	//	int num_zheng_T = num_T_total;
	//对于Potts模型最低温度一定是带有基态的

	EE[0].EE = -2 * N * N;
	EE[0].gE = log(q);
	//			num_EP[0] *= -1;
	leng_EE = 1;

	//	num_zheng_T --;

	do{
		int havefoundjg = 0;
		for (int i_T = 0; i_T < num_T_total; i_T++){
			if (num_EP[i_T]>0){
				int jg = -1;
				double P_jg = 0;
				//查看是否有overlap
				for (int i = 0; i < num_EP[i_T]; i++) {
					int j = 0;
					for (; j < leng_EE; j++) {
						if (EE[j].EE == EP[i_T][i].EE && EP[i_T][i].gE > 0) {
							break;
						}
					}
					if (j != leng_EE) {	//找到
						jg = j;
						P_jg = EP[i_T][i].gE;
						break;
					}
				}


				if (jg != -1){//存在overlap

					havefoundjg++;
					double T = get_T(T_start, T_end, i_T, num_T_total);
					double B = 1.0 / T;
					for (int i = 0; i < num_EP[i_T]; i++) {

						if (EP[i_T][i].gE <= 0){
							continue;
						}

						int j = 0;
						for (; j < leng_EE; j++){
							if (EP[i_T][i].EE == EE[j].EE)
								break;
						}
						if (j == leng_EE){
							EE[leng_EE].EE = EP[i_T][i].EE;
							EE[leng_EE].gE = EE[jg].gE + log(EP[i_T][i].gE) - log(P_jg)
								+ (-B * (EE[jg].EE - EE[leng_EE].EE));


							leng_EE++;

						}
					}
					num_EP[i_T] *= -1;
					//num_zheng_T --;
				}


			}
			//else{
			//	num_zheng_T --;
			//}

		}
		if (havefoundjg == 0){//|| num_zheng_T == 0){

			break;
		}
	} while (1);
	/*
		ofstream fout_g("g.txt");
		for(int i=0;i<leng_EE;i++){
		fout_g<<EE[i].EE<<"\t"<<EE[i].gE<<endl;
		}
		*/
	//计算g}========================================================================

	//计算F{=======================================================================
	for (int i_T = 0; i_T < num_T_total; i_T++){
		num_EP[i_T] *= -1;
		double Pmax = 0;
		double Emax = 0;
		double gE = 0;

		for (int j = 0; j < abs(num_EP[i_T]); j++){
			if (Pmax < fabs(EP[i_T][j].gE)){
				int i_E = 0;
				for (; i_E < leng_EE; i_E++){
					if (EP[i_T][j].EE == EE[i_E].EE)
						break;
				}
				if (i_E == leng_EE){
					continue;
				}
				else{
					Pmax = fabs(EP[i_T][j].gE);
					Emax = EP[i_T][j].EE;
					gE = EE[i_E].gE;
				}
			}
		}


		double T = get_T(T_start, T_end, i_T, num_T_total);
		double B = 1.0 / T;
		if (gE == 0){
			cout << "erro:no overlap E>>>T=" << T << "\ti_T=" << i_T << endl;
			continue;
		}

		if (Emax == -2 * N*N){
			double Flocal = 0;
			Flocal = (-T
				* (log(2.0 * double(m)) - log(Pmax)) + Emax)
				/ double(N * N);
			FNNf[i_T][0] += Flocal;
			FNNf[i_T][1] += Flocal*Flocal;
		}
		else{
			double Flocal = 0;
			Flocal = -T
				* (gE + (-B * Emax)
				- log(Pmax) + log(double(m)))
				/ double(N * N);
			FNNf[i_T][0] += Flocal;
			FNNf[i_T][1] += Flocal*Flocal;
			//			double F_local = FNNf[i_T][i_num_xunhuan];
			//			cout<<F_local<<endl;
		}
	}
	//计算F}=======================================================================


}

//调用overlap，计算g,计算F。利用高温信息
void P_F_high_temperature(EEx **EP, int *num_EP, int num_T_total, int Nequi, int m,
	double T_start, double T_end, double **FNNf, int xunhuan) {

	const int leng_ENE = 2 * N * N + 10;
	EEx EE[leng_ENE];
	int leng_EE = 0;
	//初始化

	for (int i = 0; i < leng_ENE; i++) {
		EE[i].EE = 1;
		EE[i].gE = 0;
	}
	leng_EE = 0;

	//step 1 overlap 将曲线连接起来
	//step 2 计算g
	//step 3 计算f


	//overlap{-----------------------------------------------------------------
	//printEP(EP,num_EP,num_T_total,1,m);

	overlap(EP, num_EP, num_T_total, m);


	//printEP(EP,num_EP,num_T_total,2,m);
	//overlap}====================================================================



	//计算g{========================================================================

	//计算无穷高温下Pmax的g{======================================================
	/*
	NOTE:
	无穷高温的Z=空间大小q^(N*N)
	g=P*Z ; log(g) = log(P)+(N*N)*log(q)
	但是P是以数量储存的，并没有归一化
	log(g) = log(P)-log(m)+(N*N)*log(q)
	*/
	//最后一个用来储存无穷温概率分布
	//找到最大概率
	for (int i = 0; i < num_EP[num_T_total - 1]; i++){
		if (EP[num_T_total - 1][i].EE == EP[num_T_total - 1][num_EP[num_T_total - 1]].EE){
			EE[0].EE = EP[num_T_total - 1][i].EE;
			EE[0].gE = log(fabs(EP[num_T_total - 1][i].gE)) - log(double(m)) + double(N*N)*log(q);
			leng_EE = 1;
			break;
		}
	}

	//Pmax g}=====================================================================

	do{
		int havefoundjg = 0;
		for (int i_T = num_T_total - 1; i_T >= 0; i_T--){
			if (num_EP[i_T]>0){
				int jg = -1;
				double P_jg = 0;
				//查看是否有overlap
				for (int i = 0; i < num_EP[i_T]; i++) {
					int j = 0;
					for (; j < leng_EE; j++) {
						if (EE[j].EE == EP[i_T][i].EE && EP[i_T][i].gE > 0) {
							break;
						}
					}
					if (j != leng_EE) {	//找到
						jg = j;
						P_jg = EP[i_T][i].gE;
						break;
					}
				}


				if (jg != -1){//存在overlap

					havefoundjg++;

					double T, B;
					if (i_T != num_T_total - 1){
						T = get_T_high_temperature(T_start, T_end, i_T, num_T_total);
						B = 1.0 / T;
					}
					else{
						T = 1.E10;
						B = 0.0;
					}
					for (int i = 0; i < num_EP[i_T]; i++) {

						if (EP[i_T][i].gE <= 0){
							continue;
						}

						int j = 0;
						for (; j < leng_EE; j++){
							if (EP[i_T][i].EE == EE[j].EE)
								break;
						}
						if (j == leng_EE){
							EE[leng_EE].EE = EP[i_T][i].EE;
							EE[leng_EE].gE = EE[jg].gE + log(EP[i_T][i].gE) - log(P_jg)
								+ (-B * (EE[jg].EE - EE[leng_EE].EE));
							leng_EE++;

						}
					}
					num_EP[i_T] *= -1;
					//num_zheng_T --;
				}


			}
			//else{
			//	num_zheng_T --;
			//}

		}
		if (havefoundjg == 0){//|| num_zheng_T == 0){

			break;
		}
	} while (1);

	//debug
	/*
	ofstream fout_g("g.txt");
	for(int i=0;i<leng_EE;i++){
	fout_g<<EE[i].EE<<"\t"<<exp(EE[i].gE)<<endl;
	}
	*/
	//debug

	//计算g}========================================================================

	//计算F{=======================================================================
	for (int i_T = 0; i_T < num_T_total - 1; i_T++){//最后一个不与要计算
		num_EP[i_T] *= -1;
		double Pmax = 0;
		double Emax = 0;
		double gE = 0;

		for (int j = 0; j < abs(num_EP[i_T]); j++){
			if (Pmax < fabs(EP[i_T][j].gE)){
				int i_E = 0;
				for (; i_E < leng_EE; i_E++){
					if (EP[i_T][j].EE == EE[i_E].EE)
						break;
				}
				if (i_E == leng_EE){
					continue;
				}
				else{
					Pmax = fabs(EP[i_T][j].gE);
					Emax = EP[i_T][j].EE;
					gE = EE[i_E].gE;
				}
			}
		}

		double T, B;
		T = get_T_high_temperature(T_start, T_end, i_T, num_T_total);
		B = 1.0 / T;

		if (gE == 0){
			cout << "erro:no overlap E>>>T=" << T << "\ti_T=" << i_T << endl;
			continue;
		}

		double Flocal = 0;
		Flocal = -T
			* (gE + (-B * Emax)
			- log(Pmax) + log(double(m)))
			/ double(N * N);
		FNNf[i_T][0] += Flocal;
		FNNf[i_T][1] += Flocal*Flocal;
		//			double F_local = FNNf[i_T][i_num_xunhuan];
		//			cout<<F_local<<endl;

	}
	//计算F}=======================================================================
}

void F_high_temperature(int m, int Nequi, double T_start, double T_end, int xunhuan, int num_T) {


	stringstream sm;
	sm << m;
	stringstream sq;
	sq << q;
	stringstream sN;
	sN << N;


	ofstream fout(("F_T_m=" + sm.str() + "_N=" + sN.str() + ".txt").c_str(), ios::app);

	double **FNNf = new double *[num_T];
#pragma omp parallel for
	for (int i = 0; i < num_T; i++) {
		FNNf[i] = new double[2];
		FNNf[i][0] = 0;
		FNNf[i][1] = 0;
	}



	time_t tstart;
	time_t tend;
	tstart = time(NULL);

	//总循环，求平均值{---------------------------------------------------------------

#pragma omp parallel for //num_threads(8)
	for (int xh = 0; xh < xunhuan; xh++) {

		printf("********************xh=%d****************\n", xh + 1);
		fflush(stdout);

		EEx **EP;
		int *num_EP = new int[num_T];
		EP = new EEx *[num_T];

		//	#pragma omp parallel for
		for (int i = 0; i < num_T; i++) {
			EP[i] = new EEx[num_E];
			num_EP[i] = 0;
			//	#pragma omp parallel for
			for (int j = 0; j < num_E; j++) {
				EP[i][j].EE = 0;
				EP[i][j].gE = 0;
			}

		}




		//采用并行，产生每个温度下的E-P数据，
		//并对能量大小进行排序，挑选出单根E-P类抛物曲线中点{-----------------------------------

		//debug
		//int num_cout = 0;
		//debug
		//#pragma omp parallel for//Attation
		for (int i = 0; i < num_T; i++) {
			//debug
			//num_cout += 1;
			//cout<< num_cout <<endl;
			//debug
			//赋值温度等参量{--------------------------------
			double T, B, P;
			int N0;
			if (i != num_T - 1){
				T = get_T_high_temperature(T_start, T_end, i, num_T);
				B = 1.0 / T;
				N0 = N0_sure(T);
				P = 1.0 - exp(-B * J);
			}
			else{
				T = get_T_high_temperature(T_start, T_end, i, num_T);
				B = 0.0;
				N0 = N0_sure(T);
				P = 0.0;
			}

			//end}==========================================

			//计算每个温度下的E-P{--------------------------------

			F_P(Nequi, m, N0, T, B, P, EP[i], num_EP[i]);

			//end}==============================================

			if (num_EP[i] == 1){
				EP[i][0].gE *= -1.0;
			}
			//printf("%d\n",num_EP[i]);
			//fflush(stdout);

		}


		//计算每个温度E-P end}=======================================================================



		//		printEP(EP,num_EP,num_T,xh,m);
		//将每个温度P-E数组传入，链接P-E图像，计算每个温度的g和F{-----------------------------------------------------------------------------


		P_F_high_temperature(EP, num_EP, num_T, Nequi, m, T_start, T_end, FNNf, xunhuan);

		//		printf("time_l=%f\n", double(time(NULL) - t1));
		//		fflush(stdout);

		//end}========================================================================



		//删除{**********************************************************************
		for (int i = 0; i < num_T; i++)
			delete[] EP[i];
		delete[]EP;

		delete[] num_EP;
		//删除}**********************************************************************
		//		printf("end--xh=%d\n",xh);

		//		fflush(stdout);
	}

	//循环结束}========================================================================


	//计算平均值和标准差{-------------------------------------------------------------
	//全温度输出

	//debug
	/*
	for (int i = 0; i < num_T-1; i++) {
	fout<<get_T_high_temperature(T_start, T_end, i, num_T)<<"\t"<<FNNf[i][0] <<endl;
	}
	*/


	for (int i = 0; i < num_T - 1; i++) {
		fout << std::setprecision(16) << FNNf[i][0] << "\t" << FNNf[i][1] << endl;
	}


	/*
	for (int i = 0; i < num_T-1; i++) {
	fout<<FNNf[i][0] <<"\t";
	}
	fout<<endl;
	for (int i = 0; i < num_T-1; i++) {
	fout<<FNNf[i][1] <<"\t";
	}
	*/

	//平均值end}=====================================================================

	for (int i = 0; i < num_T; i++)
		delete[] FNNf[i];
	delete[] FNNf;


}

/*
void F_sigma(int m, int Nequi, double T_start, double T_end, int xunhuan, int num_T,
int Sigma_i_T, double & Sigma ,double & Fvalue) {


double **FNNf = new double *[num_T];
#pragma omp parallel for
for (int i = 0; i < num_T; i++) {
FNNf[i] = new double[2];
FNNf[i][0] = 0;
FNNf[i][1] = 0;
}

//总循环，求平均值{---------------------------------------------------------------

#pragma omp parallel for num_threads(8)
for (int xh = 0; xh < xunhuan; xh++) {

printf("********************xh=%d****************\n", xh + 1);
fflush(stdout);

EEx **EP;
int *num_EP = new int[num_T];
EP = new EEx *[num_T];

//	#pragma omp parallel for
for (int i = 0; i < num_T; i++) {
EP[i] = new EEx[num_E];
num_EP[i] = 0;
//	#pragma omp parallel for
for (int j = 0; j < num_E; j++) {
EP[i][j].EE = 0;
EP[i][j].gE = 0;
}

}




//采用并行，产生每个温度下的E-P数据，
//并对能量大小进行排序，挑选出单根E-P类抛物曲线中点{-----------------------------------


#pragma omp parallel for
for (int i = 0; i < num_T; i++) {

//赋值温度等参量{--------------------------------
double T, B, P;
T = get_T(T_start, T_end, i, num_T);
B = 1.0 / T;
int N0 = N0_sure(T);
P = 1.0 - exp(-B * J);

//end}==========================================

//计算每个温度下的E-P{--------------------------------

F_P(Nequi, m, N0, T, B, P, EP[i], num_EP[i]);

//end}==============================================

}


//计算每个温度E-P end}=======================================================================




//将每个温度P-E数组传入，链接P-E图像，计算每个温度的g和F{-----------------------------------------------------------------------------

P_F(EP, num_EP, num_T, Nequi, m, T_start, T_end, FNNf, xunhuan);

//		printf("time_l=%f\n", double(time(NULL) - t1));
//		fflush(stdout);

//end}========================================================================



//删除{**********************************************************************
for(int i=0;i<num_T;i++)
delete [] EP[i];
delete []EP;

delete [] num_EP;

//删除}**********************************************************************

}

//循环结束}========================================================================

for(int i=0;i<num_T;i++){
FNNf[i][1] = sqrt(FNNf[i][1] - FNNf[i][0]*FNNf[i][0]);
}

//平均值end}=====================================================================

Sigma = FNNf[Sigma_i_T][1];
Fvalue = FNNf[Sigma_i_T][0];

for(int i=0;i<num_T;i++)
delete [] FNNf[i];
delete [] FNNf;


}
*/

void F(int m, int Nequi, double T_start, double T_end, int xunhuan, int num_T) {


	stringstream sm;
	sm << m;
	stringstream sq;
	sq << q;
	stringstream sN;
	sN << N;


	ofstream fout(("F_T_m=" + sm.str() + "_N=" + sN.str() + ".txt").c_str(), ios::app);

	double **FNNf = new double *[num_T];
#pragma omp parallel for
	for (int i = 0; i < num_T; i++) {
		FNNf[i] = new double[2];
		FNNf[i][0] = 0;
		FNNf[i][1] = 0;
	}



	time_t tstart;
	time_t tend;
	tstart = time(NULL);

	//总循环，求平均值{---------------------------------------------------------------

#pragma omp parallel for //num_threads(8)
	for (int xh = 0; xh < xunhuan; xh++) {

		printf("********************xh=%d****************\n", xh + 1);
		fflush(stdout);

		EEx **EP;
		int *num_EP = new int[num_T];
		EP = new EEx *[num_T];

		//	#pragma omp parallel for
		for (int i = 0; i < num_T; i++) {
			EP[i] = new EEx[num_E];
			num_EP[i] = 0;
			//	#pragma omp parallel for
			for (int j = 0; j < num_E; j++) {
				EP[i][j].EE = 0;
				EP[i][j].gE = 0;
			}

		}




		//采用并行，产生每个温度下的E-P数据，
		//并对能量大小进行排序，挑选出单根E-P类抛物曲线中点{-----------------------------------


		//#pragma omp parallel for
		for (int i = 0; i < num_T; i++) {

			//赋值温度等参量{--------------------------------
			double T, B, P;
			T = get_T(T_start, T_end, i, num_T);
			B = 1.0 / T;
			int N0 = N0_sure(T);
			P = 1.0 - exp(-B * J);

			//end}==========================================

			//计算每个温度下的E-P{--------------------------------

			F_P(Nequi, m, N0, T, B, P, EP[i], num_EP[i]);

			//end}==============================================

			if (num_EP[i] == 1){
				EP[i][0].gE *= -1.0;
			}
			//printf("%d\n",num_EP[i]);
			//fflush(stdout);

		}


		//计算每个温度E-P end}=======================================================================



		//		printEP(EP,num_EP,num_T,xh,m);
		//将每个温度P-E数组传入，链接P-E图像，计算每个温度的g和F{-----------------------------------------------------------------------------


		P_F(EP, num_EP, num_T, Nequi, m, T_start, T_end, FNNf, xunhuan);

		//		printf("time_l=%f\n", double(time(NULL) - t1));
		//		fflush(stdout);

		//end}========================================================================



		//删除{**********************************************************************
		for (int i = 0; i < num_T; i++)
			delete[] EP[i];
		delete[]EP;

		delete[] num_EP;
		//删除}**********************************************************************
		//		printf("end--xh=%d\n",xh);

		//		fflush(stdout);
	}

	//循环结束}========================================================================


	//计算平均值和标准差{-------------------------------------------------------------
	//全温度输出
	/*
	tend = time(NULL);
	cout << "time_total=" << (tend - tstart) << endl;
	for (int i = 0; i < num_T; i++) {
	FNNf[i][1] = sqrt(fabs(FNNf[i][1]/double(xunhuan) - FNNf[i][0]*FNNf[i][0]/double(xunhuan*xunhuan)));//将循环分为多次计算
	fout << (T_end - T_start) / double(num_T) * (i + 1.0) + T_start << "\t"
	<< FNNf[i][0]/double(xunhuan) << "\t" << FNNf[i][1] << endl;
	}
	*/
	for (int i = 0; i < num_T; i++) {
		fout << FNNf[i][0] << "\t";
	}
	fout << endl;
	for (int i = 0; i < num_T; i++) {
		fout << FNNf[i][1] << "\t";
	}
	//平均值end}=====================================================================

	for (int i = 0; i < num_T; i++)
		delete[] FNNf[i];
	delete[] FNNf;


}

//排序函数
void BubbleSort(EEx * pData, int Count) {
	double iTemp1;
	double iTemp2;
	for (int i = 0; i < Count; i++) {
		for (int j = Count - 1; j > i; j--) {
			if (pData[j].EE < pData[j - 1].EE) {
				iTemp1 = pData[j - 1].EE;
				iTemp2 = pData[j - 1].gE;
				pData[j - 1].EE = pData[j].EE;
				pData[j - 1].gE = pData[j].gE;
				pData[j].EE = iTemp1;
				pData[j].gE = iTemp2;
			}
		}
	}
}


void F_P_sigma(int Nequi, int m_p, int N0, double T, double B, double P, EEx *EP_pre_T,
	int &num_EP_prese_T, int i_m, point ** team) {



	int leng = num_EP_prese_T;

	EEx *EP_local = new EEx[num_E];
	for (int i = 0; i < num_E; i++){
		EP_local[i].gE = EP_pre_T[i].gE;
		EP_local[i].EE = EP_pre_T[i].EE;
	}

	//开始
	double HH = 0;

	if (i_m == 0){
		begin(team);
		jump(Nequi, team, P);
	}

	for (int ii = 0; ii < m_p; ii++) {

		jump(N0, team, P);

		//计算
		HH = H(team);
		int i = 0;
		for (; i < leng; i++) {
			if (EP_local[i].EE == HH)
				break;
		}

		if (i == leng) {
			EP_local[leng].EE = HH;
			EP_local[leng].gE = 1;
			leng++;
			if (leng >= num_E){
				printf("erro:num_E not enough.");
				exit(0);
			}
		}
		else {
			EP_local[i].gE++;
		}

	}



	//排序{-------------------------------------------

	BubbleSort(EP_local, leng);

	//排序end}========================================

	//找到最大能量{------------------------------------------------------

	double topE = 0;
	double topEgE = 0;

	for (int i_top = 0; i_top < leng; i_top++){
		if (topEgE < EP_local[i_top].gE){
			topE = EP_local[i_top].EE;
			topEgE = EP_local[i_top].gE;
		}
	}
	EP_local[leng].EE = topE;

	//找到最大能量}=========================================================

	//将局部变量复制到传递进来的变量{--------------------------------------------
	num_EP_prese_T = 0;

	for (int i = 0; i < leng; i++) {
		EP_pre_T[i].EE = EP_local[i].EE;
		EP_pre_T[i].gE = EP_local[i].gE;
	}
	num_EP_prese_T = leng;
	EP_pre_T[num_EP_prese_T].EE = EP_local[leng].EE;

	if (num_EP_prese_T == 0){
		printf("~~~~~~~~~~erro~~~~~~~~~~~\n");
		printf("~~~~~~~~~~~~~~~~~~~~~~~~~\n");
		fflush(stdout);
		exit(0);
	}
	//end}======================================================================

	delete[] EP_local;
}

void P_F_sigma(EEx **EP, int *num_EP, int num_T_total, int Nequi, int m,
	double T_start, double T_end, double *FNNf, int xunhuan, int i_T_sigma) {

	const int leng_ENE = 2 * N * N + 10;
	EEx EE[leng_ENE];
	int leng_EE = 0;
	//初始化

	for (int i = 0; i < leng_ENE; i++) {
		EE[i].EE = 1;
		EE[i].gE = 0;
	}
	leng_EE = 0;

	//step 1 overlap 将曲线连接起来
	//step 2 计算g
	//step 3 计算f


	//overlap{-----------------------------------------------------------------
	//	printEP(EP,num_EP,num_T_total,1);

	overlap(EP, num_EP, num_T_total, m);



	//	printEP(EP,num_EP,num_T_total,2);
	//overlap}====================================================================



	//计算g{========================================================================

	//	int num_zheng_T = num_T_total;
	//对于Potts模型最低温度一定是带有基态的

	EE[0].EE = -2 * N * N;
	EE[0].gE = log(q);
	//			num_EP[0] *= -1;
	leng_EE = 1;

	//	num_zheng_T --;

	do{
		int havefoundjg = 0;
		for (int i_T = 0; i_T < num_T_total; i_T++){
			if (num_EP[i_T]>0){
				int jg = -1;
				double P_jg = 0;
				//查看是否有overlap
				for (int i = 0; i < num_EP[i_T]; i++) {
					int j = 0;
					for (; j < leng_EE; j++) {
						if (EE[j].EE == EP[i_T][i].EE && EP[i_T][i].gE > 0) {
							break;
						}
					}
					if (j != leng_EE) {	//找到
						jg = j;
						P_jg = EP[i_T][i].gE;
						break;
					}
				}


				if (jg != -1){//存在overlap

					havefoundjg++;
					double T = get_T(T_start, T_end, i_T, num_T_total);
					double B = 1.0 / T;
					for (int i = 0; i < num_EP[i_T]; i++) {

						if (EP[i_T][i].gE <= 0){
							continue;
						}

						int j = 0;
						for (; j < leng_EE; j++){
							if (EP[i_T][i].EE == EE[j].EE)
								break;
						}
						if (j == leng_EE){
							EE[leng_EE].EE = EP[i_T][i].EE;
							EE[leng_EE].gE = EE[jg].gE + log(EP[i_T][i].gE) - log(P_jg)
								+ (-B * (EE[jg].EE - EE[leng_EE].EE));


							leng_EE++;

						}
					}
					num_EP[i_T] *= -1;
					//num_zheng_T --;
				}


			}
			//else{
			//	num_zheng_T --;
			//}

		}
		if (havefoundjg == 0){//|| num_zheng_T == 0){

			break;
		}
	} while (1);
	/*
		ofstream fout_g("g.txt");
		for(int i=0;i<leng_EE;i++){
		fout_g<<EE[i].EE<<"\t"<<EE[i].gE<<endl;
		}
		*/
	//计算g}========================================================================

	//计算F{=======================================================================

	int i_T = i_T_sigma;
	num_EP[i_T] *= -1;
	double Pmax = 0;
	double Emax = 0;
	double gE = 0;

	for (int j = 0; j < abs(num_EP[i_T]); j++){
		if (Pmax < fabs(EP[i_T][j].gE)){
			int i_E = 0;
			for (; i_E < leng_EE; i_E++){
				if (EP[i_T][j].EE == EE[i_E].EE)
					break;
			}
			if (i_E == leng_EE){
				continue;
			}
			else{
				Pmax = fabs(EP[i_T][j].gE);
				Emax = EP[i_T][j].EE;
				gE = EE[i_E].gE;
			}
		}
	}



	double T = get_T(T_start, T_end, i_T, num_T_total);
	double B = 1.0 / T;
	if (gE == 0){
		cout << "erro:no overlap E>>>T=" << T << "\ti_T=" << i_T << "\tm=" << m << endl;
	}

	if (Emax == -2 * N*N){
		double Flocal = 0;
		Flocal = (-T
			* (log(2.0 * double(m)) - log(Pmax)) + Emax)
			/ double(N * N);
		FNNf[0] += Flocal;
		FNNf[1] += Flocal*Flocal;
	}
	else{
		double Flocal = 0;
		Flocal = -T
			* (gE + (-B * Emax)
			- log(Pmax) + log(double(m)))
			/ double(N * N);
		FNNf[0] += Flocal;
		FNNf[1] += Flocal*Flocal;
		//			double F_local = FNNf[i_T][i_num_xunhuan];
		//			cout<<F_local<<endl;
	}


	//计算F}=======================================================================


}


void F_sigma(int Nequi, double T_start, double T_end, int xunhuan, int num_T,
	int Sigma_i_T, int num_m){
	//

	//由于只计算最后一个温度，所以仅仅付一个就好了
	//记得要改P_F函数
	const double Ftru = -5.20079232940099;
	double **FNNf = new double *[num_m];



	for (int i_m = 0; i_m < num_m; i_m++){
		FNNf[i_m] = new double[2];
		FNNf[i_m][0] = 0;
		FNNf[i_m][1] = 0;
	}


	time_t tstart = time(NULL);


	//总循环，求平均值{---------------------------------------------------------------
#pragma omp parallel for //num_threads(8)
	for (int xh = 0; xh < xunhuan; xh++) {
		time_t tstart_local = time(NULL);
		printf("****************xh=%d\tt=%d*************\n", xh + 1, int(tstart_local - tstart));
		fflush(stdout);

		EEx **EP;
		int *num_EP = new int[num_T];
		EP = new EEx *[num_T];

		//	#pragma omp parallel for
		for (int i = 0; i < num_T; i++) {
			EP[i] = new EEx[num_E];
			num_EP[i] = 0;
			//	#pragma omp parallel for
			for (int j = 0; j < num_E; j++) {
				EP[i][j].EE = 0;
				EP[i][j].gE = 0;
			}

		}
		int *rand_save = new int[num_T];
		point ***team = new point **[num_T];
		for (int i_team = 0; i_team < num_T; i_team++){
			rand_save[i_team] = 1;
			team[i_team] = new point *[N];
			for (int i = 0; i < N; i++) {
				team[i_team][i] = new point[N];
			}
			begin(team[i_team]);
		}
		//改变m值使得产生不同m下的概率分布{--------------------------------------------------
		for (int i_m = 0; i_m < num_m; i_m++){

			int m = int(4000 * pow(1.2, i_m));


			//每一次都得将值变正{----------------------------------------------------------
			for (int i_changabs_T = 0; i_changabs_T < num_T; i_changabs_T++){
				num_EP[i_changabs_T] = abs(num_EP[i_changabs_T]);

				for (int j_changabs = 0; j_changabs < num_E; j_changabs++){
					EP[i_changabs_T][j_changabs].gE = double(myround(fabs(EP[i_changabs_T][j_changabs].gE)));
				}
			}
			//变正end}=================================================================



			//产生不同温度下的分布{---------------------------------------------------
			//#pragma omp parallel for
			for (int i = 0; i < num_T; i++) {

				//赋值温度等参量{--------------------------------
				double T, B, P;
				T = get_T(T_start, T_end, i, num_T);
				B = 1.0 / T;
				int N0 = N0_sure(T);
				P = 1.0 - exp(-B * J);

				//end}==========================================

				//计算每个温度下的E-P{--------------------------------
				if (i_m != 0){
					seed[omp_get_thread_num()] = rand_save[i];
				}
				int m_p = i_m == 0 ? m : int(4000 * pow(1.2, i_m)) - int(4000 * pow(1.2, i_m - 1));
				point **team_local = new point *[N];
				for (int i_team = 0; i_team < N; i_team++){
					team_local[i_team] = new point[N];
					for (int j_team = 0; j_team < N; j_team++){
						team_local[i_team][j_team].cc = team[i][i_team][j_team].cc;
						team_local[i_team][j_team].x = team[i][i_team][j_team].x;
						team_local[i_team][j_team].y = team[i][i_team][j_team].y;
						team_local[i_team][j_team].num = team[i][i_team][j_team].num;

					}
				}

				F_P_sigma(Nequi, m_p, N0, T, B, P, EP[i], num_EP[i], i_m, team_local);
				rand_save[i] = seed[omp_get_thread_num()];
				for (int i_team = 0; i_team < N; i_team++){
					for (int j_team = 0; j_team < N; j_team++){
						team[i][i_team][j_team].cc = team_local[i_team][j_team].cc;
						team[i][i_team][j_team].x = team_local[i_team][j_team].x;
						team[i][i_team][j_team].y = team_local[i_team][j_team].y;
						team[i][i_team][j_team].num = team_local[i_team][j_team].num;

					}

					delete[] team_local[i_team];
				}
				delete[] team_local;
				//end}==============================================

				if (num_EP[i] == 1){
					EP[i][0].gE *= -1.0;
				}
			}
			//产生分布end}============================================================

			//		printEP(EP,num_EP,num_T,i_m,m);
			P_F_sigma(EP, num_EP, num_T, Nequi, m, T_start, T_end, FNNf[i_m], xunhuan, Sigma_i_T);

		}




		//删除{**********************************************************************
		delete[] rand_save;
		for (int i = 0; i < num_T; i++){
			delete[] EP[i];


			for (int i_N = 0; i_N < N; i_N++)
				delete[] team[i][i_N];
			delete[] team[i];

		}
		delete[]EP;

		delete[] num_EP;
		delete[] team;
		//删除}**********************************************************************


		printf("end--xh=%d\tt=%d---dt=%d---\n", xh + 1, int(time(NULL) - tstart), int(tstart_local - tstart));
		fflush(stdout);
	}
	//循环结束}========================================================================

	stringstream sdeta;
	sdeta << (T_end - T_start) / double(num_T);
	ofstream fout_Sigma(("M-Sigma_detaT=" + sdeta.str() + "_N=32_T=6.txt").c_str());
	for (int i = 0; i < num_m; i++){
		FNNf[i][1] = sqrt(fabs(FNNf[i][1] / double(xunhuan) -
			FNNf[i][0] * FNNf[i][0] / double(xunhuan*xunhuan)));//将循环分为多次计算
		fout_Sigma << int(4000 * pow(1.2, i)) << "\t" << FNNf[i][1] << "\t" << fabs(FNNf[i][0] / double(xunhuan) - Ftru) << endl;
		//	cout<<int(4000*pow(1.2,i))<<"\t"<<FNNf[i][1]<<"\t"<<FNNf[i][0]-Ftru<<endl;
	}
	fout_Sigma.close();
	//平均值end}=====================================================================


	for (int i = 0; i < num_m; i++)
		delete[] FNNf[i];
	delete[] FNNf;

}


int main(int argc, char* argv[]) //打算做q=2的potts模型的E值
{



	begin_seed(seed);

	double T, B, P;
	T = 10;   //need to be set by yourself !!!!!!!!!!!!!!!!!!!1
	B = 1.0 / T;
	int N0 = N0_sure(T);
	P = 1.0 - exp(-B * J);


	int Nequi = 1000;
	int m = 10000;

	int num_count = 0;


	point **team;
	team = new point *[N];
	for (int i = 0; i < N; i++) {
		team[i] = new point[N];
	}

	std::ofstream fout;
	fout.open("Ising_N10_T10.txt", std::ios::out);

	//开始
	double HH = 0;

	begin(team);

	jump(Nequi, team, P);


	for (int ii = 0; ii < m; ii++) {


		jump(N0, team, P);
		//if (team[0][0].num == 0){
		//	num_count += 1;
		//}


		//get configuration here which is team
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				fout << team[i][j].num << " ";
			}
		}
		fout << std::endl;


		//计算
		HH += H(team);

	}

	fout.close();
	//cout << T << " -> " << static_cast<double>(HH) / static_cast<double>(m) / static_cast<double>(N*N) << endl;

	//cout << "---------------------------" << endl;
	//cout << static_cast<double>(num_count) / static_cast<double>(m) << endl;

	//cout<<"=========================="<<endl;

	//for (int i = 0; i < N; i++) {
	//	for (int j = 0; j < N ; j++) {
	//		cout << team[i][j].num<<" ";
	//	}
	//	cout<<endl;

	//}


	//end}======================================================================
	for (int i = 0; i < N; i++)
		delete[] team[i];
	delete[] team;

	//system("qsub script_omp");
	return 0;
}