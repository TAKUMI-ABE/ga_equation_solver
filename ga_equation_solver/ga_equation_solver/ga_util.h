#ifndef GA_UTIL_H
#define GA_UTIL_H

#include<memory>

/* GAで使われるパラメータをまとめて記述するファイル */
#define F_X pow((x),2)-9 /* 探索する関数。関数への入力はx, 関数からの出力はy */
#define G_Y 1/(1+fabs(y))  /* 関数の出力値に対しての適合度。 適合度は 高いほど良いものと扱われる。f(x)=0が最良 */
#define GRAY 1             /* グレイコードかバイナリコードの指定。グレイコードなら1, バイナリコードなら0 */

/* GAのパラメータ群 */
#define MAX (5.12)        /* 扱う実数の最大値 */
#define MIN (-5.12)       /* 扱う実数の最小値 */
#define LENGTH (10)        /* 遺伝子のコード長 */
#define POP 100            /* 個体数 */
#define CODE_MAX 1         /* 各遺伝子コードの最大値。これが１ならコードは0か1になる。ビットストリングの場合は１で固定 */
#define GAP 0.9            /* 一回の生殖で子供と入れ替わる割合 */
#define ELITE_RATE 1.0     /* そのまま残る数のうち、エリートの割合 */
#define P_MUTATE 0.0      /* 突然変異率。LENGTHの逆数程度がよい */
#define P_CROSS 1.0        /* 交叉確率 */
#define GENERATION 10       /* GAを計算する世代数 */
#define SELECTION_METHOD 2 /* 1はルーレット 2はトーナメント*/
#define TOURNAMENT_SIZE 5	/* トーナメントサイズ。トーナメントの時だけ意味がある  */

/* 出力 */
#define PRINT_GROUP 1
#define PRINT_FITNESS 1

// Define the type of Genotype of GA
using gtype_t = int*;

/*
 * GAの個体を表す構造体ga_individualを宣言し,.
 * 個体は連結リストになっている.
 */
using individual_t = class ga_individual*;
//using individual_t = std::unique_ptr<struct ga_individual>;
class ga_individual {
public:
	gtype_t gtype; // Geno Type
	double ptype; // Pheno Type
	double fitness;
	individual_t next; // pointer to the next indivisual
	int rank;  // rank after sorting
	int parent1; // index of parent 1
	int parent2; // index of parent 2
	int cross_point; // crossover point
};

/* 
 * generation
 */
using ga_population_t = class ga_population*;
class ga_population {
public:
	individual_t genes;  // lead pointer to the invividual list
	double * pselect; // array of fitness
	int mutate_count;  // total number of mutation
	double max_fitness;
	double min_fitness;
	double avg_fitness;
	int population_size; // number of indivisuals
	int code_length; // length of gene
	int code_max; // maximum value of genetic Lucus
};

/*
* ptype to gtype
*/
void encode_gtype(double value, gtype_t gtype, int code_length, double min, double max)
{
	double gap = max - min;
	double remain_value = value - min; /* 値のうち、遺伝子によって表現されている部分 */
	double value_of_code;         /* その桁の遺伝子が表現する値 */
	int position = 1;
	int pre_code = 0;
	int i = 0;
	int tmp;  /* グレイコード変換用,一時保管変数 */
	while (i < code_length) {
		value_of_code = gap / pow(2, position);
		if (remain_value >= value_of_code) {
			gtype[i] = 1;
			remain_value -= value_of_code;
		}
		else {
			gtype[i] = 0;
		}
		/* グレイコードへの変換
		* バイナリコードと、元のバイナリを右に１シフトしたもののXORをとる
		*/
		if (GRAY == 1) {
			tmp = gtype[i];
			gtype[i] = (pre_code) ^ (gtype[i]);
			pre_code = tmp;
		}
		position++;
		i++;
	}
	return;
}

/*
* gtype to ptype
*/
double decode_gtype(gtype_t gtype, int code_length, double min, double max)
{
	double gap = max - min;
	double decoded_value = min;
	int position = 1;
	int pre_code = 0;                     /* １つ上位の桁のコード(バイナリ),バイナリとグレイの変換に必要 */

	int i = 0;
	/* グレイコードの解釈 */
	/* 変換されたバイナリの１つ上位の桁のコードとの排他的論理和を取る */
	if (GRAY == 1) {
		while (i < code_length) {
			pre_code = pre_code ^ gtype[i];
			if (pre_code) {
				decoded_value += gap / pow(2, position); /* 最上位から順に、最大値と最小値の差の1/2,1/4,1/8,1/16,,,となる */
			}
			position++;
			i++;
		}
	}
	/* バイナリコードの時 */
	else {
		while (i < code_length) {
			if (gtype[i]) {
				decoded_value += gap / pow(2, position); /* 最上位から順に、最大値と最小値の差の1/2,1/4,1/8,1/16,,,となる */
			}
			position++;
			i++;
		}
	}
	return decoded_value;
}

#endif // GA_UTIL_H