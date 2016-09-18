#ifndef _PARAM_H_
#define _PARAM_H_

/* GAで使われるパラメータをまとめて記述するファイル */
#define F_X pow((x),2)-10 /* 探索する関数。関数への入力はx, 関数からの出力はy */
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

#endif

