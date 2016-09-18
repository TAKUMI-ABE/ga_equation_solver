#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "ga_util.h"


	/*
	�������m�ۂ̂��߂̊֐��B�m�ۂł��Ȃ������ꍇ�̓G���[�E�I���B
	*/
	void* my_malloc(int size) {
		void* ptr = malloc(size);
		if (!ptr) {
			fprintf(stderr, "failed to maloc.\n");
			exit(1);
		}
		return ptr;
	}

	/* gtype(��`�q�R�[�h)�̎���GA�̂��߂̎��� */

	/*gtype(int�z��)�����ŏ��̃A�h���X��Ԃ�*/
	gtype_t mk_gtype(int code_length)
	{
		gtype_t gtype = (gtype_t)my_malloc(sizeof(int)*code_length);
		//gtype_t gtype = new int(sizeof(int)*code_length);
		return gtype;
	}

	/*gtype�̃��������������*/
	void free_gtype(gtype_t gtype)
	{
		free(gtype);
		return;
	}

	/*�����_����gtype�����*/
	gtype_t mk_random_gtype(int code_length, int code_max)
	{
		gtype_t ptr = mk_gtype(code_length);
		int i;
		for (i = 0; i < code_length; i++) {
			ptr[i] = rand() % (code_max + 1);
		}
		return ptr;
	}

	/* gtype �̃R�s�[*/
	void copy_gtype(gtype_t new_gtype, gtype_t old_gtype, int length) {
		int i = 0;
		for (i = 0; i < length; i++) {
			new_gtype[i] = old_gtype[i];
		}
		return;
	}

	/*gtype�̌���::�������N�������ꏊ��0����LENGTH-2��Ԃ��B
	(��)
	111111111
	000000000
	��
	100000000
	011111111
	�ƂȂ�������0��Ԃ��B
	*/
	int cross_gtype(gtype_t gtype1, gtype_t gtype2, int length) {
		int cross_point = rand() % (length - 1);
		int i = cross_point + 1;
		int tmp;
		while (i < length) {
			tmp = gtype1[i];
			gtype1[i] = gtype2[i];
			gtype2[i] = tmp;
			i++;
		}
		return cross_point;
	}

	/*gtype�̓ˑR�ψ�::�ˑR�ψق��N�������񐔂�Ԃ�
	�m��pm�œˑR�ψق��N�����B������B
	�Ⴆ��
	11101����
	10011�֕ψق�������
	3��Ԃ��B
	*/
	int  mutate_gtype(gtype_t gtype, int length, int code_max, double pm) {
		//�G���[����
		if (pm > 1 || pm < 0.0) {
			printf("%f is used for mutation probability, but this must be from 0.0 to 1.0 \n", pm);
			exit(-1);
		}
		int mutate_point = 0;
		int i = 0;
		double rm;
		for (i = 0; i < length; i++) {
			rm = (double)rand() / RAND_MAX;
			if (rm < pm) {
				gtype[i] = rand() % (code_max + 1);
				mutate_point++;
			}
		}
		return mutate_point;
	}

	/* gtype��\������ */
	void print_gtype(gtype_t gtype, int length) {
		int i = 0;
		printf("[");
		while (i < length) {
			if (gtype[i] < 10) {
				printf("%d", gtype[i]);
			}
			else {
				printf("(%d)", gtype[i]);
			}
			i++;
		}
		printf("]");
	}

	/* ���`���X�g�p�̗אڂ����v�f�̓���ւ��A�����͐擪��individual_t�̃A�h���X */
	void switch_gene(individual_t *individual) {
		individual_t tmp_ptr1 = (*individual)->next->next;
		individual_t tmp_ptr2 = (*individual)->next;
		(*individual)->next->next = *individual;
		(*individual)->next = tmp_ptr1;
		(*individual) = tmp_ptr2;
		return;
	}

	/* �̂����B�������̈�m�ہA������ */
	individual_t mk_gene(int code_length, int code_max) {
		//individual_t ptr = (individual_t)my_malloc(sizeof(struct ga_individual));
		individual_t ptr = new ga_individual;
		ptr->gtype = mk_random_gtype(code_length, code_max);
		ptr->ptype = 0;
		ptr->fitness = 0;
		ptr->next = NULL;
		ptr->parent1 = 0;
		ptr->parent2 = 0;
		ptr->cross_point = 0;
		return ptr;
	}

	/* �̂��R�s�[���� */
	void copy_gene(individual_t new_gene, individual_t old_gene, int code_length) {
		copy_gtype(new_gene->gtype, old_gene->gtype, code_length);
		new_gene->ptype = old_gene->ptype;
		new_gene->fitness = old_gene->fitness;
		new_gene->parent1 = old_gene->rank;
		new_gene->parent2 = old_gene->rank;
		new_gene->cross_point = code_length - 1;
		return;
	}

	/* �����A�ˑR�ψقŎq������� �ˑR�ψى񐔂�Ԃ� */
	int mk_children_genes(individual_t child1, individual_t child2, individual_t parent1, individual_t parent2, int code_length, int code_max, double pm) {
		int cross_point, mutateCount;
		copy_gene(child1, parent1, code_length);
		copy_gene(child2, parent2, code_length);
		child1->parent1 = parent1->rank;
		child1->parent2 = parent2->rank;
		cross_point = cross_gtype(child1->gtype, child2->gtype, code_length);
		child1->parent1 = parent1->rank;
		child1->parent2 = parent2->rank;
		child1->cross_point = cross_point;
		child2->parent1 = parent2->rank;
		child2->parent2 = parent1->rank;
		child2->cross_point = cross_point;
		mutateCount = mutate_gtype(child1->gtype, code_length, code_max, pm);
		mutateCount += mutate_gtype(child2->gtype, code_length, code_max, pm);
		return mutateCount;
	}

	/*GA�W�c�̍쐬�A���������s��*/
	ga_population_t mk_init_ga_population(int population_size, int code_length, int code_max) {
		//ga_population_t population = (ga_population_t)my_malloc(sizeof(struct ga_population));
		ga_population_t population = new struct ga_population;

		population->pselect = (double*)my_malloc(sizeof(double)*population_size);
		population->mutate_count = 0;
		population->population_size = population_size;
		population->code_length = code_length;
		population->code_max = code_max;
		individual_t list_tale;
		population->genes = mk_gene(code_length, code_max);
		list_tale = population->genes;
		int i = 1;
		for (i = 1; i < population_size; i++) {
			list_tale->next = mk_gene(code_length, code_max);
			list_tale = list_tale->next;
		}
		return population;
	}

	/* �w�肵������ch���w�肵������length�����J��Ԃ��֐� */
	/* print_population(�E)�̒��Ŏg����B */
	void print_sequence(char ch, int length) {
		int i = 0;
		for (i = 0; i < length; i++) {
			printf("%c", ch);
		}
	}

	/* �W�c��\������B */
	/* ������,���㐔,�e�̃C���f�b�N�X,�����_,gtype,ptype,fitness��\������B */
	/* �܂��A�Ō�ɓˑR�ψق̉񐔂�\������B */
	void print_population(ga_population_t population) {
		individual_t member = population->genes;
		int i = 0;
		printf("-------------------");
		print_sequence('-', LENGTH + 2);
		printf("---------------\n");
		printf("#   parents  xsite  gtype");
		print_sequence(' ', LENGTH - 3);
		printf("ptype  fitness\n");

		while (member != NULL) {
			printf("%-3d (%3d,%3d) %3d  ", i, member->parent1, member->parent2, member->cross_point);
			print_gtype(member->gtype, population->code_length);
			printf(" %+3.3f  %+3.3f\n", member->ptype, member->fitness);
			member = member->next;
			i++;
		}
		printf("total mutate %d\n", population->mutate_count);
		return;
	}

	/* �K���x���o��
	�ő�,����,�ŏ�
	CSV�`���ɂ���
	*/
	void print_fitness(ga_population_t population) {
		printf("%f, %f, %f, %f, ", population->max_fitness, population->avg_fitness, population->min_fitness, population->genes->ptype);
		print_gtype(population->genes->gtype, population->code_length);
		printf("\n");
		return;
	}

	/* GA�W�c�̌̐��`���X�ggenes�̈�l��l��fitness������
	�z��pselect�����
	1. pselect[i] = pselect[i-1]+fitness[i]
	2. pselect[i] = pselect[i]/pselect[POPULATION-1]
	*/
	void calc_pselect(ga_population_t population) {
		int i;
		population->pselect[0] = population->genes->fitness;
		individual_t gene_ptr = population->genes->next;

		for (i = 1; i < population->population_size; i++) {
			population->pselect[i] = population->pselect[i - 1] + gene_ptr->fitness;
			gene_ptr = gene_ptr->next;
		}
		for (i = 0; i < population->population_size; i++) {
			population->pselect[i] /= population->pselect[population->population_size - 1];
		}
		return;
	}

	/* �K���x�̔�r�֐� */
	/* individualA �̓K���x���Ⴏ���0�ȊO��Ԃ��A���������傫�����0��Ԃ� */
	int less_than(individual_t individualA, individual_t individualB)
	{
		return (individualA->fitness < individualB->fitness);
	}

	/* �K���x�v�Z
	gtype����ptype�ւ̕ϊ��Afitness�̌v�Z���s��
	�K���x���v�Z�����̂���K���x���ɐ��`���X�g�ɑ}������B
	param.h��F_X���֐�f(x)�Ƃ���
	param.h��G_Y��g(f(x))�Ƃ��ēK���x�Ƃ���B
	���ł�
	f(x)=0�ƂȂ�x�����߂�B
	*/
	void calc_fitness(ga_population_t population, double value_min, double value_max) {
		individual_t ptr = population->genes;
		individual_t next;
		individual_t individual_ptr = NULL;
		individual_t search_ptr = ptr;

		double x, y;
		while (ptr != NULL) {
			x = decode_gtype(ptr->gtype, population->code_length, value_min, value_max);
			ptr->ptype = x;
			y = F_X;
			ptr->fitness = G_Y;
			next = ptr->next;
			ptr->next = NULL;

			// ���`���X�g�ɓK���x���ɑ}��    
			search_ptr = individual_ptr;
			if (search_ptr == NULL || less_than(individual_ptr, ptr)) {
				ptr->next = individual_ptr;
				individual_ptr = ptr;
			}
			else {
				while (search_ptr->next != NULL) {
					if (less_than(search_ptr->next, ptr)) {
						break;
					}
					search_ptr = search_ptr->next;
				}
				ptr->next = search_ptr->next;
				search_ptr->next = ptr;
			}

			ptr = next;
		}
		population->genes = individual_ptr;
		return;
	}

	/* ���[���b�g�����ɂ��e�I�� */
	/* �e�̑I��
	�@�@1. �O����P�܂ł̗����𐶐�����B
	  �@�@2. pselect�����Đe��I��
		�@�@3. �I�΂ꂽ�e��Ԃ�
		  */
	individual_t select_parent_roulette(ga_population_t population) {
		int j = 0;
		double r;
		individual_t parent;
		r = (double)rand() / RAND_MAX;
		parent = population->genes;

		/* pselect[j-1]<r<pselect[j]�̂Ƃ��Aj�Ԗڂ��e�ɂȂ� */
		while (r > population->pselect[j]) {
			parent = parent->next;
			j++;
		}
		return parent;
	}

	/* �g�[�i�����g�����ɂ��e�I�� */
	individual_t select_parent_tournament(ga_population_t population, int tournament_size) {
		int pop = population->population_size;
		int i, j, r, min = pop;
		individual_t min_selected = NULL;
		individual_t ptr;


		for (i = 0; i < tournament_size; i++) {
			r = rand() % pop;
			if (min > r)
			{
				min = r;
			}
		}
		ptr = population->genes;
		for (j = 0; j < min; j++) {
			ptr = ptr->next;
		}
		min_selected = ptr;
		return min_selected;
	}

	/* �e�̂̑I���Aparam.h��SELECTION_METHOD�ɂ����
	���[���b�g�I�����g�[�i�����g�I�����s�� */
	individual_t select_parent(ga_population_t population) {
		individual_t parent;
		switch (SELECTION_METHOD) {
		case 1:
			parent = select_parent_roulette(population);
			break;
		case 2:
			parent = select_parent_tournament(population, TOURNAMENT_SIZE);
			break;
		default:
			fprintf(stderr, "invalid number on SELECTION_METHOD\n");
			exit(1);
		}
		return parent;
	}

	/*
	�K���x���ɕ��񂾐��`���X�g����
	�ő�l�A�ŏ��l�A���ϒl���L�^�A���ԕt��
	*/
	void normalize_population(ga_population_t population) {
		int i;
		individual_t tmp;

		tmp = population->genes;
		population->max_fitness = population->genes->fitness; /* �擪�̓K���x���ő�K���x */
		double avg = 0.0;
		/* ���ԕt��*/
		for (i = 0; i < population->population_size; i++) {
			avg += tmp->fitness;
			tmp->rank = i;
			population->min_fitness = tmp->fitness;
			tmp = tmp->next;
		}
		avg = avg / population->population_size;
		population->avg_fitness = avg;

		return;
	}

	/* �V��������̐���
	new_population�̃������̈�͂��łɊm�ۂ��Ă���Ƃ���
	�K���\�[�g�ς݂�population��n������
	*/
	void generate_population(ga_population_t new_population, ga_population_t old_population, double gap, double elite_rate, double mutate_prob, double crossover_prob) {
		int num_of_remain = (int)(old_population->population_size*(1 - gap)); /*�e���ォ��R�s�[���鐔*/
		int num_of_elite = (int)(num_of_remain*elite_rate);  /*�R�s�[�g�̂����G���[�g�̐�*/
		int generated;
		double rand_double;
		individual_t old_gene = old_population->genes;
		individual_t new_gene = new_population->genes;

		/* �e�I���e�[�u�������� */
		calc_pselect(old_population);

		/* �G���[�g�헪 �e����ł̏�ʈ�萔�͂��̂܂܎q���ɂȂ� */
		for (generated = 0; generated < num_of_elite; generated++) {
			copy_gene(new_gene, old_gene, old_population->code_length);
			old_gene = old_gene->next;
			new_gene = new_gene->next;
		}
		/* �G���[�g�ȊO�̂��̂܂܎q���ɂȂ�g */
		for (; generated < num_of_remain; generated++) {
			copy_gene(new_gene, select_parent(old_population), old_population->code_length);
			new_gene = new_gene->next;
		}

		new_population->mutate_count = 0;
		/* �����E�ˑR�ψق�K������g */
		/* �c��̐�����̎��́A������ˑR�ψقŎq������� */
		if ((old_population->population_size - generated) % 2 == 1) {
			copy_gene(new_gene, select_parent(old_population), old_population->code_length);
			new_population->mutate_count += mutate_gtype(new_gene->gtype, old_population->code_length, old_population->code_max, mutate_prob);
			new_gene = new_gene->next;
			generated++;
		}

		/* �����E�ˑR�ψق����� */
		for (; generated < old_population->population_size; generated += 2) {
			rand_double = (double)rand() / RAND_MAX;
			/* ��������Ƃ� */
			if (rand_double < crossover_prob) {
				new_population->mutate_count += mk_children_genes(new_gene, new_gene->next, select_parent(old_population), select_parent(old_population), old_population->code_length, old_population->code_max, mutate_prob);
				new_gene = new_gene->next->next;
			}
			/* �������Ȃ��Ƃ� */
			else {
				copy_gene(new_gene, select_parent(old_population), old_population->code_length);
				new_population->mutate_count += mutate_gtype(new_gene->gtype, old_population->code_length, old_population->code_max, mutate_prob);
				new_gene = new_gene->next;
				copy_gene(new_gene, select_parent(old_population), old_population->code_length);
				new_population->mutate_count += mutate_gtype(new_gene->gtype, old_population->code_length, old_population->code_max, mutate_prob);
				new_gene = new_gene->next;
			}
		}
		return;
	}


	/* main�֐� */
	/* GA�̎��s */
	int main() {
		/* �����Ɉ�����^���� */
		srand(time(NULL));

		ga_population_t parent_group = mk_init_ga_population(POP, LENGTH, CODE_MAX);
		ga_population_t child_group = mk_init_ga_population(POP, LENGTH, CODE_MAX);
		int i;
		if (PRINT_FITNESS == 1) printf("#generation,max_fitness, avg_fitness, min_fitness, best_individual_ptype, best_individual_gtype\n");

		for (i = 0; i <= GENERATION; i++) {
			// �W�c�̓K���x���v�Z���A���`���X�g�����
			calc_fitness(parent_group, MIN, MAX);
			// �ő�l�E�ŏ��l�A
			normalize_population(parent_group);

			// ���݂̐���̕\��
			if (PRINT_GROUP == 1) {
				print_population(parent_group);
			}
			if (PRINT_FITNESS == 1) {
				printf("%3d, ", i);
				print_fitness(parent_group);
			}

			// ���݂̐���parent_group���玟����child_group�����B
			generate_population(child_group, parent_group, GAP, ELITE_RATE, P_MUTATE, P_CROSS);

			// ��������ւ���B
			parent_group = child_group;
		}

		return 0;
	}