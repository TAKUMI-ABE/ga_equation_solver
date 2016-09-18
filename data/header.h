#ifndef _HEADER_H_
#define _HEADER_H_

/* intの配列をGAの遺伝子型として別名を付ける。 */
/*   プログラム中ではgtype_tを使う。*/
typedef int* gtype_t;

/* GAの個体を表す構造体ga_individualを宣言し、*/
/* そのポインタ型をga_individual_tとする。 */
/* 個体は連結リストになっている。 */
typedef struct ga_individual* individual_t;
struct ga_individual{
  gtype_t gtype;  /*遺伝子型 int配列の先頭のポインタ*/
  double ptype;  /* 表現型 */
  double fitness; /* 適応度*/
  individual_t next; /* 連結リストでの次の個体 */
  int rank;  /* 連結リストのなかでの順位(ソート後) */
  int parent1; /* 交叉での親1のインデックス */
  int parent2; /* 交叉での親2のインデックス */
  int cross_point; /* 交叉したポイント */
};

/* 集団を表す構造体を宣言し、*/
/*   そのポインタ型をga_population_tとする。 */
typedef struct ga_population* ga_population_t;
struct ga_population {
  individual_t genes;  /* 個体の連結リスト先頭のポインタ */
  double * pselect; /* 適合度の配列 */
  int mutate_count;  /* 突然変異回数の合計 */
  double max_fitness;  /* 適合度の最大値 */
  double min_fitness;  /* 適合度の最小値 */
  double avg_fitness;  /* 適合度の平均値 */
  int population_size; /* 集団の個体数 */
  int code_length; /* 遺伝子長 */
  int code_max; /* 各遺伝子座の最大値,ビットストリングの場合は1 */
};

#endif
