#ifndef GA_UTIL_H
#define GA_UTIL_H

#include<memory>

/* GA�Ŏg����p�����[�^���܂Ƃ߂ċL�q����t�@�C�� */
#define F_X pow((x),2)-9 /* �T������֐��B�֐��ւ̓��͂�x, �֐�����̏o�͂�y */
#define G_Y 1/(1+fabs(y))  /* �֐��̏o�͒l�ɑ΂��Ă̓K���x�B �K���x�� �����قǗǂ����̂ƈ�����Bf(x)=0���ŗ� */
#define GRAY 1             /* �O���C�R�[�h���o�C�i���R�[�h�̎w��B�O���C�R�[�h�Ȃ�1, �o�C�i���R�[�h�Ȃ�0 */

/* GA�̃p�����[�^�Q */
#define MAX (5.12)        /* ���������̍ő�l */
#define MIN (-5.12)       /* ���������̍ŏ��l */
#define LENGTH (10)        /* ��`�q�̃R�[�h�� */
#define POP 100            /* �̐� */
#define CODE_MAX 1         /* �e��`�q�R�[�h�̍ő�l�B���ꂪ�P�Ȃ�R�[�h��0��1�ɂȂ�B�r�b�g�X�g�����O�̏ꍇ�͂P�ŌŒ� */
#define GAP 0.9            /* ���̐��B�Ŏq���Ɠ���ւ�銄�� */
#define ELITE_RATE 1.0     /* ���̂܂܎c�鐔�̂����A�G���[�g�̊��� */
#define P_MUTATE 0.0      /* �ˑR�ψٗ��BLENGTH�̋t�����x���悢 */
#define P_CROSS 1.0        /* �����m�� */
#define GENERATION 10       /* GA���v�Z���鐢�㐔 */
#define SELECTION_METHOD 2 /* 1�̓��[���b�g 2�̓g�[�i�����g*/
#define TOURNAMENT_SIZE 5	/* �g�[�i�����g�T�C�Y�B�g�[�i�����g�̎������Ӗ�������  */

/* �o�� */
#define PRINT_GROUP 1
#define PRINT_FITNESS 1

// Define the type of Genotype of GA
using gtype_t = int*;

/*
 * GA�̌̂�\���\����ga_individual��錾��,.
 * �̂͘A�����X�g�ɂȂ��Ă���.
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
	double remain_value = value - min; /* �l�̂����A��`�q�ɂ���ĕ\������Ă��镔�� */
	double value_of_code;         /* ���̌��̈�`�q���\������l */
	int position = 1;
	int pre_code = 0;
	int i = 0;
	int tmp;  /* �O���C�R�[�h�ϊ��p,�ꎞ�ۊǕϐ� */
	while (i < code_length) {
		value_of_code = gap / pow(2, position);
		if (remain_value >= value_of_code) {
			gtype[i] = 1;
			remain_value -= value_of_code;
		}
		else {
			gtype[i] = 0;
		}
		/* �O���C�R�[�h�ւ̕ϊ�
		* �o�C�i���R�[�h�ƁA���̃o�C�i�����E�ɂP�V�t�g�������̂�XOR���Ƃ�
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
	int pre_code = 0;                     /* �P��ʂ̌��̃R�[�h(�o�C�i��),�o�C�i���ƃO���C�̕ϊ��ɕK�v */

	int i = 0;
	/* �O���C�R�[�h�̉��� */
	/* �ϊ����ꂽ�o�C�i���̂P��ʂ̌��̃R�[�h�Ƃ̔r���I�_���a����� */
	if (GRAY == 1) {
		while (i < code_length) {
			pre_code = pre_code ^ gtype[i];
			if (pre_code) {
				decoded_value += gap / pow(2, position); /* �ŏ�ʂ��珇�ɁA�ő�l�ƍŏ��l�̍���1/2,1/4,1/8,1/16,,,�ƂȂ� */
			}
			position++;
			i++;
		}
	}
	/* �o�C�i���R�[�h�̎� */
	else {
		while (i < code_length) {
			if (gtype[i]) {
				decoded_value += gap / pow(2, position); /* �ŏ�ʂ��珇�ɁA�ő�l�ƍŏ��l�̍���1/2,1/4,1/8,1/16,,,�ƂȂ� */
			}
			position++;
			i++;
		}
	}
	return decoded_value;
}

#endif // GA_UTIL_H