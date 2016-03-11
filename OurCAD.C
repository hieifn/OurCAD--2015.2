/*
Elementos aceitos e linhas do netlist:

Resistor:  R<nome> <no+> <no-> <resistencia>
VCCS:      G<nome> <io+> <io-> <vi+> <vi-> <transcondutancia>
VCVC:      E<nome> <vo+> <vo-> <vi+> <vi-> <ganho de tensao>
CCCS:      F<nome> <io+> <io-> <ii+> <ii-> <ganho de corrente>
CCVS:      H<nome> <vo+> <vo-> <ii+> <ii-> <transresistencia>
Fonte I:   I<nome> <io+> <io-> <corrente>
Fonte V:   V<nome> <vo+> <vo-> <tensao> <parâmetros>
Amp. op.:  O<nome> <vo1> <vo2> <vi1> <vi2>
Indutor:   L<nome> <no+> <no-> <indutancia> [Corrente Inicial]
Capacitor: C<nome> <no+> <no-> <Capacitância> [Tensão Inicial]

Parâmetros para as fontes de tensão (V):
    Fonte DC : DC <valor>
    Fonte senoidal: SIN <nível contínuo> <amplitude> <frequencia> <atraso> <amortecimento> <defasagem>
    Fonte pulso: PULSE <amplitude1> <amplitude2> <atraso> <tempo de subida> <tempo de descida> <tempo ligada> <período> <número de cíclos>

As fontes F e H tem o ramo de entrada em curto
O amplificador operacional ideal tem a saida suspensa
Os nos podem ser nomes
*/

#define versao "1.0j - 13/02/2016"
#include <stdio.h>
#include <conio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <time.h>
#include <math.h>
#define MAX_LINHA 80
#define MAX_NOME 11
#define MAX_ELEM 50
#define MAX_NOS 50
#define TOLG 1e-9
#define TOLG2 3e-20
#define MAX_ERRO_NR 1e-29
//#define DEBUG
#define MAX_STEPS 4

typedef struct elemento { /* Elemento do netlist */
  char nome[MAX_NOME],type[MAX_NOME];
  double valor,ini, param1, param2, param3, param4, param5, param6, param7,g1,g2,g3,i1,i2,i3;
  int a,b,c,d,x,y;

} elemento;

elemento netlist[MAX_ELEM]; /* Netlist */

int
  tryAgain=1,
  RaphsonCount=0,
  repete=1,
  goToNewton,
  useInicialConditions = 2,
  quant,
  intSteps,
  order,
  save = 3, /* variavel para armazenamento  */
  icc = 0, /* Correntes capacitor */
  ns,
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
  timeA=0,
  iStepSize= TOLG2,
  finalTime,
  stepSize,
  g,
  z,
  pulseOffTime,
  pulseRealTime,
  Yc[MAX_STEPS+2][MAX_NOS+1], /* Essa matriz ira armazenar os valores das correntes nos capacitores */
  Ys[MAX_STEPS+2][MAX_NOS+1], /* Essa matriz irá armazenar até no máximo 3 valores passados da analise */
  Ynr[MAX_NOS+1][MAX_NOS+2],
  NRCompare[MAX_NOS+1],
  nrErro[MAX_NOS+1],
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

void AdamsMoltonL (int i, int save) /* ADMO do Indutor Completo! */
{
    if (save == 3){
        if(useInicialConditions==1){netlist[i].ini = 0;}
        z=(netlist[i].ini*(netlist[i].valor/iStepSize));
        g=netlist[i].valor/iStepSize;
    }
    else if (order == 1){ /* PERFEITO */
        z=((netlist[i].valor/stepSize)*(Ys[save+1][netlist[i].x]));
        g=(netlist[i].valor/stepSize);
    }
    else if (order == 2){/* PERFEITO */
        z=((((2.0*netlist[i].valor)/stepSize)*(Ys[save+1][netlist[i].x]))+(Ys[save+1][netlist[i].a]-Ys[save+1][netlist[i].b]));
        g=((2.0*netlist[i].valor)/stepSize);
     //   printf("z: %g g: %g\n ",z,g);
     //   getch();
    }
    else if (order == 3){ /* PERFEITO */
        z=( ((12.0/5.0)*((Ys[save+1][netlist[i].x])*((netlist[i].valor)/(stepSize)))) + ((8.0/5.0)*(Ys[save+1][netlist[i].a]-Ys[save+1][netlist[i].b])) - ((1.0/5.0)*(Ys[save+2][netlist[i].a]-Ys[save+2][netlist[i].b])));
        g=((12.0/5.0)*((netlist[i].valor)/stepSize));
      //  printf("z: %g g: %g",z,g);
      //  getch();
    }
    else if (order == 4){ /* Perfeito */
        z=((((Ys[save+1][netlist[i].x])*((netlist[i].valor)/stepSize))*(24.0000/9.0000))+((19.0000/9.0000)*(Ys[save+1][netlist[i].a]-Ys[save+1][netlist[i].b]))-((5.0000/9.0000)*(Ys[save+2][netlist[i].a]-Ys[save+2][netlist[i].b]))+((1.0000/9.0000)*(Ys[save+3][netlist[i].a]-Ys[save+3][netlist[i].b])));
        g=((24.0000/9.0000)*((netlist[i].valor)/stepSize));
     //   printf("z: %g g: %g\n ",z,g);
     //   getch();
    }
}

void AdamsMoltonC (int i, int save)
{
    if (save == 3){
        if(useInicialConditions==1){netlist[i].ini = 0;}   /* iStpeSize - StepSize inicial para começar a analise */
        z=(netlist[i].ini*(netlist[i].valor/iStepSize)); /* iStpeSize - StepSize inicial para começar a analise */
        g=((netlist[i].valor)/(iStepSize));

    }
    else if (order == 1){ /* PERFEITO */
        z=((Ys[save+1][netlist[i].a]-Ys[save+1][netlist[i].b])*(netlist[i].valor/stepSize));
        g=((netlist[i].valor)/(stepSize));
    }
    else if (order == 2){ /* PERFEITO */
        z=(((Ys[save+1][netlist[i].a]-Ys[save+1][netlist[i].b])*((2.0*netlist[i].valor)/stepSize))+Yc[save+1][netlist[i].x]);
        g=(2.0*(netlist[i].valor)/(stepSize));
    //    printf("z :%g g: %g\n",z,g);
    //    getch();
    }
    else if (order == 3){ /* PERFEITO */
        z=(((Ys[save+1][netlist[i].a]-Ys[save+1][netlist[i].b])*(((12.0/5.0)*netlist[i].valor)/stepSize))+ ((8.0/5.0)*Yc[save+1][netlist[i].x]) - ((1.0/5.0)*Yc[save+2][netlist[i].x])  );
        g=((12.0/5.0)*((netlist[i].valor)/(stepSize)));
//        printf("save0: %g save2: %g save3: %g\n",(Yc[save+1][netlist[i].x]),Yc[save+2][netlist[i].x],Yc[save+3][netlist[i].x]);
//        getch();
    }
    else if (order == 4){ /* ERRO NA ESTAMPA, PROVAVELMENTE SINAL */
        z=(  ((Ys[save+1][netlist[i].a]-Ys[save+1][netlist[i].b])*((24.0000/9.0000)*((netlist[i].valor)/(stepSize)))) + ((19.0000/9.0000)*Yc[save+1][netlist[i].x]) - ((5.0000/9.0000)*Yc[save+2][netlist[i].x]) + ((1.0000/9.0000)*Yc[save+3][netlist[i].x]));
        g=((24.0000/9.0000)*((netlist[i].valor)/(stepSize)));
 //       printf("save0: %g save2: %g save3: %g\n",(Yc[save+1][netlist[i].x]),Yc[save+2][netlist[i].x],Yc[save+3][netlist[i].x]);
 //       getch();
    }
}

int main(void)
{
  system("cls");
  srand( (unsigned)time(NULL) );
  printf("Programa de analise no tempo, pelo metodo de Adams-Molton\n");
  printf("Desenvolvido por : Igor F. Nascimento, Eduardo Naslausky e Alan Carpilovsky\n"); /*Galera troquem deoois !*/
  printf("Versao %d\n",rand());
  //printf("Versao %s\n",versao);
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
    printf("%lg %lg %s %i %s %i\n",finalTime,stepSize,method,intSteps,uic,quant);/* Debug - Igor */
    order=atoi(method+4); /* Tem que ser o 4 pq o 5 é o endOfString ADMO"N"  */
    ne--;
    }
    else if (tipo=='$'){
      sscanf(p,"%10s%10s%10s%10s%lg%lg%lg",na,nb,nc,nd,
                                            &netlist[ne].param1,
                                            &netlist[ne].param2,
                                            &netlist[ne].param3);
        printf("%s %s controle1:%s controle2:%s %g\n",netlist[ne].nome,na,nb,netlist[ne].valor);
        netlist[ne].a=numero(na);
        netlist[ne].b=numero(nb);
        netlist[ne].c=numero(nc);
        netlist[ne].d=numero(nd);

    }
    else if (tipo=='N'){
      sscanf(p,"%10s%10s%lg%lg%lg%lg%lg%lg%lg%lg",na,nb,&netlist[ne].param1,
                                            &netlist[ne].param2,
                                            &netlist[ne].param3,
                                            &netlist[ne].param4,
                                            &netlist[ne].param5,
                                            &netlist[ne].param6,
                                            &netlist[ne].param7,
                                            &netlist[ne].valor);

      netlist[ne].i1 =((((netlist[ne].param4-netlist[ne].param2)/(netlist[ne].param3-netlist[ne].param1))*-netlist[ne].param1)+netlist[ne].param2);
      netlist[ne].g1 =((netlist[ne].param4-netlist[ne].param2)/(netlist[ne].param3-netlist[ne].param1));
      netlist[ne].i2 =((((netlist[ne].param6-netlist[ne].param4)/(netlist[ne].param5-netlist[ne].param3))*-netlist[ne].param3)+netlist[ne].param4);
      netlist[ne].g2 =((netlist[ne].param6-netlist[ne].param4)/(netlist[ne].param5-netlist[ne].param3));
      netlist[ne].i3 =((((netlist[ne].valor-netlist[ne].param6)/(netlist[ne].param7-netlist[ne].param5))*-netlist[ne].param5)+netlist[ne].param6);
      netlist[ne].g3 =((netlist[ne].valor-netlist[ne].param6)/(netlist[ne].param7-netlist[ne].param5));

      netlist[ne].a=numero(na);
      netlist[ne].b=numero(nb);
      printf("%s %s %s P1<%g,%g> P2<%g,%g> P3<%g,%g> P4<%g,%g>\n ",netlist[ne].nome,na,nb,netlist[ne].param1,netlist[ne].param2,netlist[ne].param3,netlist[ne].param4,netlist[ne].param5,netlist[ne].param6,netlist[ne].param7,netlist[ne].valor);
    //  printf("i1  %g  g1 %g i2 %g g2 %g i3 %g g3 %g \n", netlist[ne].i1, netlist[ne].g1, netlist[ne].i2, netlist[ne].g2, netlist[ne].i3, netlist[ne].g3);
    } /*
     Conferir depois se tá lendo certinho tb*/


    else if (tipo=='L' || tipo=='C'){
      if ((quant=sscanf(p,"%10s%10s%lg IC=%lg" ,na,nb,&netlist[ne].valor, &netlist[ne].ini))!=4){
        netlist[ne].ini=0; /* caso UIC não seja especificada */
      };
      printf("%s %s %s %g %g quant =%i\n",netlist[ne].nome,na,nb,netlist[ne].valor,netlist[ne].ini,quant);
      netlist[ne].a=numero(na);
      netlist[ne].b=numero(nb);
    }
    else if (tipo=='V') {
      sscanf(p,"%10s%10s%10s%lg%lg%lg%lg%lg%lg%lg%lg",na,nb,netlist[ne].type,
                                            &netlist[ne].valor, /* Amplitude 1*/
                                            &netlist[ne].param1,/* Amplitude2*/
                                            &netlist[ne].param2, /* Atraso*/
                                            &netlist[ne].param3, /* Tempo de subida*/
                                            &netlist[ne].param4, /* Tempo de descida*/
                                            &netlist[ne].param5, /* Tempo ligada*/
                                            &netlist[ne].param6, /* Periodo */
                                            &netlist[ne].param7);/* n de ciclos*/
      printf("%s %s %s %g\n",netlist[ne].nome,na,nb,netlist[ne].valor);
   // printf("%s %g %g %g %g %g %g %g \n", netlist[ne].type, netlist[ne].param1,netlist[ne].param2,netlist[ne].param3,netlist[ne].param4,netlist[ne].param5,netlist[ne].param6,netlist[ne].param7);
      netlist[ne].a=numero(na);
      netlist[ne].b=numero(nb);
    }
    else if (tipo=='R' || tipo=='I') {
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
    else if (tipo=='C'){ /* Salva no .x onde estarão salvas as correntes dos capacitores. */
      icc++;
      netlist[i].x=icc;
    }
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
    else if (tipo=='N') {
      printf("%s %d %d %g %g %g %g %g %g\n",netlist[i].nome,netlist[i].a,netlist[i].b,netlist[i].i1,netlist[i].g1,netlist[i].i2,netlist[i].g2,netlist[i].i3,netlist[i].g3);
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
	printf("Leio e copio");
  getch();
  return 0;
}

