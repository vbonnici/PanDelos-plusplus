#ifndef PANDELOS_PLUSPLUS_CUSTOM_LIST_H
#define PANDELOS_PLUSPLUS_CUSTOM_LIST_H

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

struct elemento {
    unsigned int valore;
    unsigned int valore2;
    struct elemento *successivo;
};

struct elemento *crea_lista() {

    struct elemento *p, *punt;

    p = (struct elemento *) malloc(sizeof(struct elemento));

    punt = p;

    punt->successivo = NULL; // marcatore fine lista

    return p;
}

void stampa_lista(struct elemento *testa) {
    struct elemento *e;

    for (e = testa; e != NULL; e = e->successivo) {
        printf("%d, %d \n", e->valore, e->valore2);
    }

    printf("\n");
}

/* ***************START: PILA *******************/

void metti_su_pila(struct elemento **cima, unsigned int valore, unsigned int valore2) {
    struct elemento *nuovo;

    nuovo = (struct elemento *) malloc(sizeof(struct elemento));
    nuovo->valore = valore;
    nuovo->valore2 = valore2;
    nuovo->successivo = *cima;      // sposto la cima in basso

    *cima = nuovo;                  // in cima c'è il nuovo elemento
}

#endif //PANDELOS_PLUSPLUS_CUSTOM_LIST_H
