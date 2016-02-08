/*
Elementos aceitos e linhas do netlist:

Resistor:  R<nome> <no+> <no-> <resistencia>
VCCS:      G<nome> <io+> <io-> <vi+> <vi-> <transcondutancia>
VCVC:      E<nome> <vo+> <vo-> <vi+> <vi-> <ganho de tensao>
CCCS:      F<nome> <io+> <io-> <ii+> <ii-> <ganho de corrente>
CCVS:      H<nome> <vo+> <vo-> <ii+> <ii-> <transresistencia>
Fonte I:   I<nome> <io+> <io-> <corrente>
Fonte V:   V<nome> <vo+> <vo-> <tensao>
Amp. op.:  O<nome> <vo1> <vo2> <vi1> <vi2>
Indutor:   L<nome> <no+> <no-> <indutancia> [Corrente Inicial]
Capacitor: C<nome> <no+> <no-> <Capacitância> [Tensão Inicial]

As fontes F e H tem o ramo de entrada em curto
O amplificador operacional ideal tem a saida suspensa
Os nos podem ser nomes
*/

#define versao "1.0j - 26/11/2015"
#include <stdio.h>
#include <conio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#define MAX_LINHA 80
#define MAX_NOME 11
#define MAX_ELEM 50
#define MAX_NOS 50
#define TOLG 1e-9
#define TOLG2 3e-2
#define DEBUG
#define MAX_STEPS 4

typedef struct elemento { /* Elemento do netlist */
  char nome[MAX_NOME];
  double valor,ini;
  int a,b,c,d,x,y;
} elemento;

elemento netlist[MAX_ELEM]; /* Netlist */

int
  useInicialConditions = 2,
  quant,
  intSteps,
  order,
  ne, /* Elementos */
  nv, /* Variaveis */
  nn, /* Nos */
  i,j,k;

char
/* Foram colocados limites nos formatos de leitura para alguma protecao
   contra excesso de caracteres nestas variaveis */
  nomearquivo[MAX_LINHA+1],
  tipo,
  na[MAX_NOME],nb[MAX_NOME],nc[MAX_NOME],nd[MAX_NOME],
  lista[MAX_NOS+1][MAX_NOME+2], /*Tem que caber jx antes do nome */
  txt[MAX_LINHA+1],
  method[6],
  uic[10],
  *p;
FILE *arquivo;

double
  iStepSize= TOLG2,
  finalTime,
  stepSize,
  g,
  z,
  Ys[MAX_NOS+1][MAX_STEPS+1], /* Essa matriz irá armazenar até no máximo 5 valores passados da analise, mas necessitamos salvar tudo em um arquivo */
  Yn[MAX_NOS+1][MAX_NOS+2];

/* Resolucao de sistema de equacoes lineares.
   Metodo de Gauss-Jordan com condensacao pivotal */
int resolversistema(void)
{
  int i,j,l, a;
  double t, p;

  for (i=1; i<=nv; i++) {
    t=0.0;
    a=i;
    for (l=i; l<=nv; l++) {
      if (fabs(Yn[l][i])>fabs(t)) {
	a=l;
	t=Yn[l][i];
      }
    }
    if (i!=a) {
      for (l=1; l<=nv+1; l++) {
	p=Yn[i][l];
	Yn[i][l]=Yn[a][l];
	Yn[a][l]=p;
      }
    }
    if (fabs(t)<TOLG) {
      printf("Sistema singular\n");
      return 1;
    }
    for (j=nv+1; j>0; j--) {  /* Basta j>i em vez de j>0 */
      Yn[i][j]/= t;
      p=Yn[i][j];
      if (p!=0)  /* Evita operacoes com zero */
        for (l=1; l<=nv; l++) {
	  if (l!=i)
	    Yn[l][j]-=Yn[l][i]*p;
        }
    }
  }
  return 0;
}

/* Rotina que conta os nos e atribui numeros a eles */
int numero(char *nome)
{
  int i,achou;

  i=0; achou=0;
  while (!achou && i<=nv)
    if (!(achou=!strcmp(nome,lista[i]))) i++;
  if (!achou) {
    if (nv==MAX_NOS) {
      printf("O programa so aceita ate %d nos\n",nv);
      exit(1);
    }
    nv++;
    strcpy(lista[nv],nome);
    return nv; /* novo no */
  }
  else {
    return i; /* no ja conhecido */
  }
}

int main(void)
{
  system("cls");
  printf("Programa de analise no tempo, pelo metodo de Adams-Molton\n");
  printf("Por Transões da UFRJ\n"); /*Galera troquem deoois !*/
  printf("Versao %s\n",versao);
 denovo:
  /* Leitura do netlist */
  ne=0; nv=0; strcpy(lista[0],"0");
  printf("Nome do arquivo com o netlist (ex: mna.net): ");
  scanf("%50s",nomearquivo);
  arquivo=fopen(nomearquivo,"r");
  if (arquivo==0) {
    printf("Arquivo %s inexistente\n",nomearquivo);
    goto denovo;
  }
  printf("Lendo netlist:\n");
  fgets(txt,MAX_LINHA,arquivo);
  printf("Titulo: %s",txt);
  while (fgets(txt,MAX_LINHA,arquivo)) {
    ne++; /* Nao usa o netlist[0] */
    if (ne>MAX_ELEM) {
      printf("O programa so aceita ate %d elementos\n",MAX_ELEM);
      exit(1);
    }
    txt[0]=toupper(txt[0]);
    tipo=txt[0];
    sscanf(txt,"%10s",netlist[ne].nome);
    p=txt+strlen(netlist[ne].nome); /* Inicio dos parametros */
    /* O que e lido depende do tipo */
    if (tipo=='.'){
      if ((quant=sscanf(p,"%lg%lg%6s%i%4s",&finalTime,&stepSize,method,&intSteps,uic))!=5){
        useInicialConditions = 1; /* 1 = Não usar Condições Iniciais ; 2 = Usar Condições Iniciais */
      }
    
    #ifdef DEBUG  
    printf("%lg %lg %s %i %s quant=%i\n",finalTime,stepSize,method,intSteps,uic,quant);     /* Debug - Igor ;)  */
    #endif
    
    order=atoi(method+4); /* Tem que ser o 4 pq o 5 é o endOfString ADMO"N"  */
    ne--;
    }
    else if (tipo=='L' || tipo=='C'){
      if ((quant=sscanf(p,"%10s%10s%lg IC=%lg" ,na,nb,&netlist[ne].valor, &netlist[ne].ini))!=4){
        netlist[ne].ini=0; /* caso UIC não seja especificada */
      };
      printf("%s %s %s %g %g quant =%i\n",netlist[ne].nome,na,nb,netlist[ne].valor,netlist[ne].ini,quant);
      netlist[ne].a=numero(na);
      netlist[ne].b=numero(nb);
    }
    else if (tipo=='R' || tipo=='I' || tipo=='V') {
      sscanf(p,"%10s%10s%lg",na,nb,&netlist[ne].valor);
      printf("%s %s %s %g\n",netlist[ne].nome,na,nb,netlist[ne].valor);
      netlist[ne].a=numero(na);
      netlist[ne].b=numero(nb);
    }
    else if (tipo=='G' || tipo=='E' || tipo=='F' || tipo=='H') {
      sscanf(p,"%10s%10s%10s%10s%lg",na,nb,nc,nd,&netlist[ne].valor);
      printf("%s %s %s %s %s %g\n",netlist[ne].nome,na,nb,nc,nd,netlist[ne].valor);
      netlist[ne].a=numero(na);
      netlist[ne].b=numero(nb);
      netlist[ne].c=numero(nc);
      netlist[ne].d=numero(nd);
    }
    else if (tipo=='O') {
      sscanf(p,"%10s%10s%10s%10s",na,nb,nc,nd);
      printf("%s %s %s %s %s\n",netlist[ne].nome,na,nb,nc,nd);
      netlist[ne].a=numero(na);
      netlist[ne].b=numero(nb);
      netlist[ne].c=numero(nc);
      netlist[ne].d=numero(nd);
    }
    else if (tipo=='*') { /* Comentario comeca com "*" */
      printf("Comentario: %s",txt);
      ne--;
    }
    else {
      printf("Elemento desconhecido: %s\n",txt);
      getch();
      exit(1);
    }
  }
  fclose(arquivo);
  /* Acrescenta variaveis de corrente acima dos nos, anotando no netlist */
  nn=nv;
  for (i=1; i<=ne; i++) {
    tipo=netlist[i].nome[0];
    if (tipo=='V' || tipo=='E' || tipo=='F' || tipo=='O' || tipo=='L') {
      nv++;
      if (nv>MAX_NOS) {
        printf("As correntes extra excederam o numero de variaveis permitido (%d)\n",MAX_NOS);
        exit(1);
      }
      strcpy(lista[nv],"j"); /* Tem espaco para mais dois caracteres */
      strcat(lista[nv],netlist[i].nome);
      netlist[i].x=nv;
    }
    /* Retirei o Else if que estava só para o L e juntei com o acima,
    dos outros elementos que tem variáveis de corrente adicionadas
    de forma igual  */
    
    else if (tipo=='H') {
      nv=nv+2;
      if (nv>MAX_NOS) {
        printf("As correntes extra excederam o numero de variaveis permitido (%d)\n",MAX_NOS);
        exit(1);
      }
      strcpy(lista[nv-1],"jx"); strcat(lista[nv-1],netlist[i].nome);
      netlist[i].x=nv-1;
      strcpy(lista[nv],"jy"); strcat(lista[nv],netlist[i].nome);
      netlist[i].y=nv;
    }
  }
  getch();
  /* Lista tudo */
  printf("Variaveis internas: \n");
  for (i=0; i<=nv; i++)
    printf("%d -> %s\n",i,lista[i]);
  getch();
  printf("Netlist interno final\n");
  for (i=1; i<=ne; i++) {
    tipo=netlist[i].nome[0];
    if (tipo=='R' || tipo=='I' || tipo=='V') {
      printf("%s %d %d %g\n",netlist[i].nome,netlist[i].a,netlist[i].b,netlist[i].valor);
    }
    else if (tipo=='C' || tipo=='L'){
      printf("%s %d %d %g IC = %g\n",netlist[i].nome,netlist[i].a,netlist[i].b,netlist[i].valor,netlist[i].ini);
    }
    else if (tipo=='G' || tipo=='E' || tipo=='F' || tipo=='H') {
      printf("%s %d %d %d %d %g\n",netlist[i].nome,netlist[i].a,netlist[i].b,netlist[i].c,netlist[i].d,netlist[i].valor);
    }
    else if (tipo=='O') {
      printf("%s %d %d %d %d\n",netlist[i].nome,netlist[i].a,netlist[i].b,netlist[i].c,netlist[i].d);
    }
    if (tipo=='V' || tipo=='E' || tipo=='F' || tipo=='O'){
      printf("Corrente jx: %d\n",netlist[i].x);
      }
    else if (tipo=='L'){
      printf("Corrente %s: %d\n",lista[netlist[i].x],netlist[i].x);
    }
    else if (tipo=='H'){
      printf("Correntes %s: %d e %s: %d\n",lista[netlist[i].x],netlist[i].x,lista[netlist[i].y],netlist[i].y);
    }
  }
  getch();
  /* Monta o sistema nodal modificado */
  printf("O circuito tem %d nos, %d variaveis e %d elementos\n",nn,nv,ne);
  getch();
  /* Zera sistema */
  for (i=0; i<=nv; i++) {
    for (j=0; j<=nv+1; j++)
      Yn[i][j]=0;
  }
  /* Monta estampas */ /* ANTES DE TUDO FAZER ANALISE DE ORDEM 1 COM STEP MENOR QUE O NORMAL PARA ACHAR O V(t0) */
  for (i=1; i<=ne; i++) {
    tipo=netlist[i].nome[0];
    if(tipo=='L'){                          /* Para ajudar :   j(t0),v(t0) - Fonte de corrente fixa que vai de A -> B */                                                                             /* (∆t/nL(...)) - Resistor (INVERTIDO) 1/R */
      if(useInicialConditions==1){netlist[i].ini = 0;}
        z=netlist[i].ini; /* Atenção aqui; X denota as condições iniciais. Caso não tenha foi SETADO como 0 */
        g=iStepSize/netlist[i].valor; /* iStpeSize - StepSize inicial para começar a analise */  /* 1 - j(t0+∆t) ≈ j(t0) + (∆t/L)v(t0 + ∆t) */
        Yn[netlist[i].a][netlist[i].x]+=1;                                   /* 2 - j(t0+∆t) ≈ j(t0) + (∆t/2L)(v(t0) + v(t0 + ∆t)) */
        Yn[netlist[i].b][netlist[i].x]-=1;                                   /* 3 - j(t0+∆t) ≈ j(t0) + (∆t/12L)(8v(t0) - v(t0 - ∆t) + 5v(t0 + ∆t)) */
        Yn[netlist[i].x][netlist[i].a]-=1;                                   /* 4 - j(t0+∆t) ≈ j(t0) + (∆t/24L)(19v(t0) - 5v(t0 - ∆t) + 3v(t0 + ∆t)+ v(t0 - 2∆t)) */
        Yn[netlist[i].x][netlist[i].b]+=1;
        Yn[netlist[i].x][netlist[i].x]+=g;
        Yn[netlist[i].x][nv+1]+=z;          /* Fonte de corrente sendo Adicionada */
    }
    else if(tipo=='C'){                          /* Para ajudar :   j(t0),v(t0) - Fonte de corrente fixa que vai de A -> B */                                                                             /* (∆t/nL(...)) - Resistor (INVERTIDO) 1/R */
      if(useInicialConditions==1){netlist[i].ini = 0;}   /* iStpeSize - StepSize inicial para começar a analise */
        z=(netlist[i].ini*(netlist[i].valor/iStepSize)); /* iStpeSize - StepSize inicial para começar a analise */
        g=(netlist[i].valor/iStepSize);
        Yn[netlist[i].a][netlist[i].a]+=g;
        Yn[netlist[i].b][netlist[i].b]+=g;
        Yn[netlist[i].a][netlist[i].b]-=g;
        Yn[netlist[i].b][netlist[i].a]-=g;
        Yn[netlist[i].a][nv+1]+=z; /* DETALHE PARA O SINAL Fonte de corrente sendo Adicionada */
        Yn[netlist[i].b][nv+1]-=z; /* DETALHE PARA O SINAL Fonte de corrente sendo Adicionada */
    }
    else if (tipo=='R') {
      g=1/netlist[i].valor;
      Yn[netlist[i].a][netlist[i].a]+=g;
      Yn[netlist[i].b][netlist[i].b]+=g;
      Yn[netlist[i].a][netlist[i].b]-=g;
      Yn[netlist[i].b][netlist[i].a]-=g;
    }
    else if (tipo=='G') {
      g=netlist[i].valor;
      Yn[netlist[i].a][netlist[i].c]+=g;
      Yn[netlist[i].b][netlist[i].d]+=g;
      Yn[netlist[i].a][netlist[i].d]-=g;
      Yn[netlist[i].b][netlist[i].c]-=g;
    }
    else if (tipo=='I') {
      g=netlist[i].valor;
      Yn[netlist[i].a][nv+1]-=g;
      Yn[netlist[i].b][nv+1]+=g;
    }
    else if (tipo=='V') {
      Yn[netlist[i].a][netlist[i].x]+=1;
      Yn[netlist[i].b][netlist[i].x]-=1;
      Yn[netlist[i].x][netlist[i].a]-=1;
      Yn[netlist[i].x][netlist[i].b]+=1;
      Yn[netlist[i].x][nv+1]-=netlist[i].valor;
    }
    else if (tipo=='E') {
      g=netlist[i].valor;
      Yn[netlist[i].a][netlist[i].x]+=1;
      Yn[netlist[i].b][netlist[i].x]-=1;
      Yn[netlist[i].x][netlist[i].a]-=1;
      Yn[netlist[i].x][netlist[i].b]+=1;
      Yn[netlist[i].x][netlist[i].c]+=g;
      Yn[netlist[i].x][netlist[i].d]-=g;
    }
    else if (tipo=='F') {
      g=netlist[i].valor;
      Yn[netlist[i].a][netlist[i].x]+=g;
      Yn[netlist[i].b][netlist[i].x]-=g;
      Yn[netlist[i].c][netlist[i].x]+=1;
      Yn[netlist[i].d][netlist[i].x]-=1;
      Yn[netlist[i].x][netlist[i].c]-=1;
      Yn[netlist[i].x][netlist[i].d]+=1;
    }
    else if (tipo=='H') {
      g=netlist[i].valor;
      Yn[netlist[i].a][netlist[i].y]+=1;
      Yn[netlist[i].b][netlist[i].y]-=1;
      Yn[netlist[i].c][netlist[i].x]+=1;
      Yn[netlist[i].d][netlist[i].x]-=1;
      Yn[netlist[i].y][netlist[i].a]-=1;
      Yn[netlist[i].y][netlist[i].b]+=1;
      Yn[netlist[i].x][netlist[i].c]-=1;
      Yn[netlist[i].x][netlist[i].d]+=1;
      Yn[netlist[i].y][netlist[i].x]+=g;
    }
    else if (tipo=='O') {
      Yn[netlist[i].a][netlist[i].x]+=1;
      Yn[netlist[i].b][netlist[i].x]-=1;
      Yn[netlist[i].x][netlist[i].c]+=1;
      Yn[netlist[i].x][netlist[i].d]-=1;
    }
#ifdef DEBUG
    /* Opcional: Mostra o sistema apos a montagem da estampa */
    printf("Sistema apos a estampa de %s\n",netlist[i].nome);
    for (k=1; k<=nv; k++) {
      for (j=1; j<=nv+1; j++)
        if (Yn[k][j]!=0) printf("%+3.1f ",Yn[k][j]);
        else printf(" ... ");
      printf("\n");
    }
    getch();
#endif
  }
  /* Resolve o sistema */
  if (resolversistema()) {
    getch();
    exit;
  }
#ifdef DEBUG
  /* Opcional: Mostra o sistema resolvido */
  printf("Sistema resolvido:\n");
  for (i=1; i<=nv; i++) {
      for (j=1; j<=nv+1; j++)
        if (Yn[i][j]!=0) printf("%+3.1f ",Yn[i][j]);
        else printf(" ... ");
      printf("\n");
    }
  getch();
#endif
  /* Mostra solucao */
  printf("Solucao:\n");
  strcpy(txt,"Tensao");
  for (i=1; i<=nv; i++) {
    if (i==nn+1) strcpy(txt,"Corrente");
    printf("%s %s: %g\n",txt,lista[i],Yn[i][nv+1]);
  }
  getch();
  return 0;
}

