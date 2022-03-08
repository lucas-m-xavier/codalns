#include "pctf.h"

#include <conio.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include <stdio.h>
#include <string.h>

//#define DBG // habilita o modo DEBUG (exibe na tela as melhoras na FO)
#define RND // habilita múltiplas (NUM_EXE) execuções com sementes aleatórias para a instância

#define MAX(x,y) (((x) > (y)) ? (x) : (y))
#define MIN(x,y) (((x) < (y)) ? (x) : (y))


//============================== DADOS DE ENTRADA ==============================
static char INST[100]           = "C:\\Users\\L\\Documents\\UFES\\ALNS_PROJETO\\Instancias\\ES.txt";  // arquivo com a instância (estado)
static int ALFA                 = 30;       // limite de contadores instalados (% do total de arestas)
static int BETA                 = 30;       // limite de faixas cobertas (% do total de faixas)
static int NUM_EXE              = 5;        // número de execuções do método
static double MAX_TIME          = 40;      // tempo máximo de execução (segundos)
static int MAX_ITERACOES        = 100000;   // número maximo de iterações (ignnorado pelo tempo)
static int SEG                  = 25;       // tamanho do segmento--
static double FAT_REACAO        = 0.05;     // fator de reação
static int S1                   = 50;       // score se a melhor solução é superada
static int S2                   = 10;       // score se a solução corrente é superada
static int S3                   = 7;        // score se uma solução pior é aceita
static double TEMP_INICIAL      = 10000;    // temperatura inicial--- 10 100 1000
static double TEMP_CONGELAMENTO = 0.01;     // temperatura de congelamento
static double TAXA_RESFRIAMENTO = 0.95;    // taxa de resfriamento
static double WST_PMT           = 3;        // parâmetro da heurística "worst" 2 3 
static double SHW_PMT           = 3;        // parâmetro da heurística "shaw" 2 3
static double PER_REM           = 0.1;     // percentual máximo de contadores removidos 10 - 60
//==============================================================================


//================================== PRINCIPAL =================================
int main(int argc, char *argv[]) { 
    FILE *f;
    FILE* fil = fopen("calibracao_alns.txt", "a");
    Solucao sol;
    char arq[100],dir[100];
    int bstSol;
    double solMed,medSolAva,temMed,desvio;
    #ifdef RND
        srand(time(NULL));
    #else
        NUM_EXE = 1;
    #endif	
    //-----------------------
    // parâmetros
    if (argc > 1) {
        strcpy(INST, argv[1]);
        ALFA = 30;
        BETA = 10;
        NUM_EXE = 1;
        TEMP_INICIAL = atoi(argv[2]);
        TEMP_CONGELAMENTO = 0.01; 
        TAXA_RESFRIAMENTO = atof(argv[3]);
        MAX_ITERACOES = atoi(argv[4]);
        MAX_TIME = 10.0;
        SEG = atoi(argv[5]);
        FAT_REACAO = atoi(argv[6]);
        S1 = atoi(argv[7]);
        S2 = atoi(argv[8]);
        S3 = atoi(argv[9]);
    }
    //-----------------------
    // dados de entrada
    sprintf(arq,"%s",INST);
    strcpy(dir,"Solucoes");
    lerInstancia(arq);
    montarRede();
    // executar o ALNS
    bstSol = 0;
    solMed = medSolAva = temMed = desvio = 0.0;
        for (int r = 1; r <= NUM_EXE; r++) {
            //printf("\n\n>>> Resolvendo a instancia %s ALFA = %d e BETA = %d - rodada %d\n\n", INST, ALFA, BETA, r);
            execALNS(sol);
            #ifdef DBG
            sprintf(arq, "%s\\sol-%s-%da-%db-%d.txt", dir, INST, ALFA, BETA, r);
            escreverResultado(sol, arq);
            #endif
            if (sol.numParCob > bstSol) {
                #ifdef DBG
                sprintf(arq, "%s\\bstSol-%s-%da-%db.txt", dir, INST, ALFA, BETA);
                escreverResultado(sol, arq);
                #endif
                bstSol = sol.numParCob;
            }
            
            solMed += sol.numParCob;
            medSolAva += solAva_;
            temMed += bstTime_;
            #ifdef DBG
            f = fopen("saida-full.txt", "at");
            fprintf(f, "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t\t%d\t%d\t%d\t%.2f\t\t%.2f\n",
            INST, ALFA, BETA, numNos_, numAre_, numTotAre_, numTotFai_, numPar_, solAva_, sol.numConIns, sol.numFaiCob, sol.numParCob, bstTime_, excTime_);
            fclose(f);
            #endif
        }
    //-----------------------
    // calcular as médias
    solMed = solMed / (double)NUM_EXE;
    medSolAva = medSolAva / (double)NUM_EXE;
    temMed = temMed / (double)NUM_EXE;
    desvio = ((bstSol - solMed) / bstSol) * 100;
    fprintf(fil, "==========================================\nFO: %.2f\ttempo: %f\n", solMed, temMed);
    fprintf(fil, "max_iterações: %d\tseg: %d\tqtd_pior: %d\tqtd_shaw: %d\n==========================================\n", MAX_ITERACOES, SEG);
    fclose(fil);
    #ifdef DBG
    f = fopen("saida.txt", "at");
    fprintf(f, "%s\t%d\t%d\t%d\t%.2f\t%.2f\t%.2f\t%.2f\n", INST, ALFA, BETA, bstSol, solMed, desvio, temMed, medSolAva);
    fclose(f);
    #endif
    //-----------------------
    // limpar memória
    limparMemoria();
    printf("%d", bstSol);

    return 0;
}
//==============================================================================

//================================== SOLUÇÃO ===================================

//------------------------------------------------------------------------------
void criarSolucao(Solucao &s) {
    int foAntes,bstFO,bstPos;
    int vetAre[MAX_ARE];
    memset(vetAre,0,sizeof(vetAre));
    s.numConIns = s.numFaiCob = 0;
    s.numAreCom = numAre_;
    for(int i = 0; i < numAre_; i++) {
       vetAre[i] = 1;
       s.numConIns += vetArestas_[i].nCon;
       s.numFaiCob += vetArestas_[i].nFai;
    }
    s.numParCob = maxPar_;
    while((s.numConIns > maxContReal_) || (s.numFaiCob > maxFaix_)) {
        bstFO = bstPos = -1;

        for(int i = 0; i < numAre_; i++)
         if(vetAre[i] == 1) {
            foAntes = s.numParCob;
 		    vetAre[i] = 0;
            calcParCob2(s, vetAre);

            if(s.numParCob > bstFO) {
                bstFO = s.numParCob;
                bstPos = i;
            }
            vetAre[i] = 1;
            s.numParCob = foAntes;
         } 
       vetAre[bstPos] = 0;
       s.numConIns -= vetArestas_[bstPos].nCon;
       s.numFaiCob -= vetArestas_[bstPos].nFai;
       s.numAreCom--;
       s.numParCob = bstFO;
    }

    int aux = 0, aux2 = 0;
    for (int i = 0; i < numAre_; i++) {
        if (vetAre[i] == 1) {
            s.vetIDCon[aux] = i;
            aux++;
        }
        else {
            s.vetIDSem[aux2] = i;
            aux2++;
        }
    }
    s.numAreSem = numAre_ - s.numAreCom;
    /*
    printf("\n%d\n", numAre_);
    for (int cont = 0; cont < numAre_; cont++) {
        printf("%d ", vetAre[cont]);
    }
    printf("\n%d\n", s.numAreCom);
    for (int cont = 0; cont < s.numAreCom; cont++) {
        printf("%d ", s.vetIDCon[cont]);
    }
    printf("\n%d\n", s.numAreSem);
    for (int cont = 0; cont < s.numAreSem; cont++) {
        printf("%d ", s.vetIDSem[cont]);
    }
    */
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
void calcParCob(Solucao &s)
{
 int res = 0;
 int vetAre[MAX_ARE];
 memset(vetAre, 0, sizeof(vetAre));
 for (int cont = 0; cont < s.numAreCom; cont++) {
     vetAre[s.vetIDCon[cont]] = 1;
 }
 #pragma omp parallel for reduction(+:res)
 for(int n = 0; n < numPar_-1; n++)
   res += calcPar(s,n, vetAre);
 s.numParCob = maxPar_ - res;
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

void calcParCob2(Solucao& s, const int *vetAre)
{
    int res = 0;
#pragma omp parallel for reduction(+:res)
    for (int n = 0; n < numPar_ - 1; n++)
        res += calcPar(s, n, vetAre);
    s.numParCob = maxPar_ - res;
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
int calcPar(Solucao &s,const int &n1, const int *vetAre)
{
 int qtd;
 int vetBFSVis[MAX_NOS]; 
 int vetBFSFila[MAX_NOS];
 int res,u,v,iniFila,tamFila;
 res = 0;
 qtd = numPar_ - (n1+1);
 memset(vetBFSVis,0,sizeof(vetBFSVis));
 vetBFSVis[n1] = 1;
 vetBFSFila[0] = n1;
 iniFila = 0;
 tamFila = 1;
 while(tamFila > 0)
  {
   u = vetBFSFila[iniFila];
   iniFila++;
   tamFila--;
   for(int i = 0; i < vetQtdAdj_[u]; i++)
     if(vetAre[matAreAdj_[u][i]] == 0)
      {
       v = matNosAdj_[u][i];
       if(!vetBFSVis[v])
        {
         if((v > n1) && (v < numPar_))
          {
           res++;
           qtd--;
           if(qtd == 0)
             return res;
          }
         vetBFSVis[v] = 1;
	        vetBFSFila[iniFila+tamFila] = v; // circular: vetBFSFila_[(iniFila+tamFila)%MAX_NOS] = v; - não foi necessário
         tamFila++;
        }
      }
  }
	return res;
}
//------------------------------------------------------------------------------

//==============================================================================


//============================== ENTRADA E SAÍDA ===============================

//------------------------------------------------------------------------------
void lerInstancia(char *arq)
{
 bool flag;
 int qtdNos,no1,no2;
 FILE *f = fopen(arq,"r");
 fscanf(f,"%d %d %d %d %d %d\n",&numNos_,&numAre_,&numPar_,&numAreVir_,&numTotAre_,&numTotFai_);
 vetIdPar_ = new int[numPar_];
 vetIdNos_ = new int[numNos_];
 for(int i = 0; i < numPar_; i++)
  {
   fscanf(f,"%d\n",&vetIdPar_[i]);
   vetIdNos_[i] = vetIdPar_[i];
  }
 qtdNos = numPar_;
 vetArestas_ = new Aresta[numAre_];
 for(int i = 0; i < numAre_; i++)
  {
   fscanf(f,"%d %d %d %d %d\n",&vetArestas_[i].id,&no1,&no2,&vetArestas_[i].nCon,&vetArestas_[i].nFai);
   flag = 1;
   for(int j = 0; j < qtdNos; j++)
     if(no1 == vetIdNos_[j])
      {
       vetArestas_[i].no1 = j;
       flag = 0;
       break;
      }
   if(flag)
    {
     vetIdNos_[qtdNos] = no1;
     vetArestas_[i].no1 = qtdNos;  
     qtdNos++;
    }
   flag = 1;
   for(int j = 0; j < qtdNos; j++)
     if(no2 == vetIdNos_[j])
      {
       vetArestas_[i].no2 = j;
       flag = 0;
       break;
      }
   if(flag)
    {
     vetIdNos_[qtdNos] = no2;
     vetArestas_[i].no2 = qtdNos;  
     qtdNos++;
    }
  }
 fclose(f);
 maxPar_ = (numPar_*(numPar_-1))/2; 
	maxCont_ = floor(numTotAre_*ALFA/100.0);
	maxFaix_ = floor(numTotFai_*BETA/100.0);
 maxContReal_ = MIN(maxCont_,maxFaix_/2);
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
void escreverSolucao(Solucao &s,FILE *f)
{
 if(f == NULL)
   f = stdout;
	fprintf(f,"\n< ---------------------------- SOLUCAO ----------------------------- >\n");
	fprintf(f,"Numero de postos de contagem instalados.............: %d (max = %d)\n",s.numConIns,maxCont_);
	fprintf(f,"Numero de faixas cobertas...........................: %d (max = %d)\n",s.numFaiCob,maxFaix_);
	fprintf(f,"Numero de pares OD cobertos (FO)....................: %d\n",s.numParCob);
 fprintf(f,"\nVetor solucao: [ ");
 for(int i = 0; i < numAre_; i++)
   //fprintf(f,"%d ",s.vetAre[i]);
 fprintf(f,"]\n\n");
 fprintf(f,"-----------------------\n>>> ARESTAS PARA INSTALACAO DOS POSTOS DE CONTAGEM:\n");
 //for(int i = 0; i < numAre_; i++)
   //if(s.vetAre[i] == 1)
 	   //fprintf(f,"%d (%d)\n",vetArestas_[i].id,vetArestas_[i].nCon); 
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
void escreverResultado(Solucao &s,char *path)
{
	FILE *f = fopen(path,"w");	
	fprintf(f,"< ---------------------------- PROBLEMA ---------------------------- >\n");
	fprintf(f,"Instancia...........................................: %s\n",INST);
    fprintf(f,"Reduzida............................................: NAO\n");
	fprintf(f,"Limite de contadores instalados.....................: %d (%d%% arestas sem reducao)\n",maxCont_,ALFA);
	fprintf(f,"Limite de faixas cobertas...........................: %d (%d%% faixas sem reducao)\n",maxFaix_,BETA);
	fprintf(f,"Numero de nos.......................................: %d\n",numNos_);
	fprintf(f,"Numero de arestas...................................: %d\n",numAre_);
	fprintf(f,"Numero de arestas virtuais..........................: %d\n",numAreVir_);
	fprintf(f,"Numero total de arestas (sem redução)...............: %d\n",numTotAre_);
	fprintf(f,"Numero total de faixas (sem redução)................: %d\n",numTotFai_);
	fprintf(f,"Numero de nos que definem os pares..................: %d\n",numPar_);
	fprintf(f,"Numero total de pares...............................: %d\n",maxPar_);
	fprintf(f,"\n< --------------------------- PARAMETROS --------------------------- >\n");
	fprintf(f,"Temperatura inicial.................................: %.2f\n",TEMP_INICIAL);
	fprintf(f,"Temperatura de congelamento.........................: %.3f\n",TEMP_CONGELAMENTO);
	fprintf(f,"Taxa de resfriamento................................: %.3f\n",TAXA_RESFRIAMENTO);
	fprintf(f,"N. maximo de iteracoes..............................: %d\n", MAX_ITERACOES);
	fprintf(f,"Tempo maximo de processamento.......................: %.2f segundos!\n",MAX_TIME);
	fprintf(f,"\n< --------------------------- RESULTADOS --------------------------- >\n");
	fprintf(f,"Custo da solucao inicial............................: %d\n",foIni_);
	fprintf(f,"Custo da solucao final..............................: %d\n",foFin_);
	fprintf(f,"Numero de solucoes avaliadas........................: %d\n",solAva_);
	fprintf(f,"Tempo para encontrar a melhor solucao (BST).........: %.2f segundos!\n",bstTime_);
	fprintf(f,"Tempo de execucao do metodo.........................: %.2f segundos!\n",excTime_);
 escreverSolucao(s,f);
 fclose(f);
}
//------------------------------------------------------------------------------

//==============================================================================


//================================= AUXILIARES =================================

//------------------------------------------------------------------------------
void montarRede()
{
 int aux;
 int vetNos[MAX_NOS];
 int vetAre[MAX_NOS];
 matNosAdj_ = new int*[numNos_];
 matAreAdj_ = new int*[numNos_];
 vetQtdAdj_ = new int[numNos_];
 for(int i = 0; i < numNos_; i++)
  {
   vetQtdAdj_[i] = 0;
   for(int j = 0; j < numAre_;j++)
    {
     aux = -1;
     if(vetArestas_[j].no1 == i)
       aux = vetArestas_[j].no2;
     else if(vetArestas_[j].no2 == i)
       aux = vetArestas_[j].no1;
     if(aux != -1)
      {
       vetNos[vetQtdAdj_[i]] = aux;
       vetAre[vetQtdAdj_[i]] = j;
       vetQtdAdj_[i]++;
      }
    }
   matNosAdj_[i] = new int[vetQtdAdj_[i]];
   matAreAdj_[i] = new int[vetQtdAdj_[i]];
   for(int j = 0; j < vetQtdAdj_[i];j++)
    {
     matNosAdj_[i][j] = vetNos[j];
     matAreAdj_[i][j] = vetAre[j];
    }
  }
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
void limparMemoria()
{
 delete[] vetIdPar_;
 delete[] vetIdNos_;
 delete[] vetArestas_;
 delete[] vetQtdAdj_;
 for(int i = 0; i < numNos_; i++)
   delete[] matNosAdj_[i];
 delete[] matNosAdj_;
 for(int i = 0; i < numNos_; i++)
   delete[] matAreAdj_[i];
 delete[] matAreAdj_;
}
//------------------------------------------------------------------------------

//==============================================================================

//================================= ALNS =================================

//------------------------------------------------------------------------------

void calculaCusto(Solucao& s) {
    memset(cost, 0, sizeof(cost));
    int foAntes = s.numParCob;

    for (int i = 0; i < s.numAreCom; i++) {
        removeElemento(s, i);
        calcParCob(s);
        cost[i] = foAntes - s.numParCob;
        adicionaElemento(s, s.numAreSem - 1);
        s.numParCob = foAntes;
    }
}

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

void ordenaCusto(Solucao& s) {
    int i, j, aux1, aux2;

    for (i = 1; i < s.numAreCom; i++)
        for (j = s.numAreCom - 1; j >= i; j--)
            if(cost[j] > cost[j - 1]) {
                aux1 = cost[j];
                aux2 = s.vetIDCon[j];
                cost[j] = cost[j - 1];
                s.vetIDCon[j] = s.vetIDCon[j - 1];
                cost[j- 1] = aux1;
                s.vetIDCon[j - 1] = aux2;
            }
}

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
void ordenaCustoShaw(int aux, Solucao& s) {
    int i, j, aux1, aux2;
    for (i = 1; i < s.numAreCom; i++)
        for (j = s.numAreCom - 1; j >= i; j--)
            if(abs(aux - cost[j]) < abs(aux - cost[j - 1])) {
                aux1 = cost[j];
                aux2 = s.vetIDCon[j];
                cost[j] = cost[j - 1];
                s.vetIDCon[j] = s.vetIDCon[j - 1];
                cost[j- 1] = aux1;
                s.vetIDCon[j - 1] = aux2;
            }
}

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

void removeElemento(Solucao& s, int aux) {
    s.vetIDSem[s.numAreSem] = s.vetIDCon[aux];
    s.numAreSem++;
    s.numConIns -= vetArestas_[s.vetIDCon[aux]].nCon;
    s.numFaiCob -= vetArestas_[s.vetIDCon[aux]].nFai;
    s.numAreCom--;
    s.vetIDCon[aux] = s.vetIDCon[s.numAreCom];
}

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

void adicionaElemento(Solucao& s, int aux) {
    s.vetIDCon[s.numAreCom] = s.vetIDSem[aux];
    s.numAreCom++;
    s.numConIns += vetArestas_[s.vetIDSem[aux]].nCon;
    s.numFaiCob += vetArestas_[s.vetIDSem[aux]].nFai;
    s.numAreSem--;
    s.vetIDSem[aux] = s.vetIDSem[s.numAreSem];
}

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//Testar em instancias maiores (calculo de cujsto dentro do for)
void piorRemocao(Solucao& s, int aux) {
    int i;

    calculaCusto(s);
    ordenaCusto(s);

    for (int cont = 0; cont < aux; cont++) {
        i = floor(s.numAreCom * (pow((rand() / (RAND_MAX + 1.0f)), WST_PMT)));
        removeElemento(s, i);
    }
}

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

void aleRemocao(Solucao& s, int aux) {
    int val;
    calculaCusto(s);
    for (int cont = 0; cont < aux; cont++) {
        val = rand() % s.numAreCom;
        removeElemento(s, val);
    }
}

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

bool estaEmauxL(int n, int lista[]) {
    for(int cont = 0; cont < numAre_; cont++){
        if(lista[cont] == n) return true;
    }

    return false;
}

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

void shawRemocao(Solucao& s, int aux) {
    int auxL[MAX_ARE], r = rand() % s.numAreCom, i;

    removeElemento(s, r);
    auxL[0] = r;

    for (int cont = 1; cont < aux; cont++) {
        calculaCusto(s);
        r = rand() % cont;

        ordenaCustoShaw(cost[auxL[r]], s);

        i = floor(numAre_ * (pow((rand() / (RAND_MAX + 1.0f)), SHW_PMT)));

        if (estaEmauxL(i, auxL)) removeElemento(s, pos[i + 1]);
        else removeElemento(s, pos[i]);

        auxL[cont] = pos[i];
    }
}

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

void reconstruirSolucaoAle(Solucao& s) {
    int aux;
    while ((s.numConIns < maxContReal_) && (s.numFaiCob < maxFaix_)) {
        aux = rand() % s.numAreSem;
        adicionaElemento(s, aux);
    }
    calcParCob(s);
}

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

void reconstruirSolucaoPar(Solucao& s) {
    double x;
    int foAntes = s.numParCob, bstFO, bstPosA, bstPosB;

    while ((s.numConIns < maxContReal_) && (s.numFaiCob < maxFaix_)) {
        bstFO = bstPosA = bstPosB = -1;
        for (int cont = 0; cont < s.numAreSem - 1; cont++) {
            for (int cont1 = cont + 1; cont1 < s.numAreSem; cont1++) {
                adicionaElemento(s, cont);
                adicionaElemento(s, cont1);
                calcParCob(s);

                if (s.numParCob > bstFO) {
                    bstFO = s.numParCob;
                    bstPosA = cont;
                    bstPosB = cont1;
                }

                removeElemento(s, s.numAreCom - 1);
                removeElemento(s, s.numAreCom - 1);
            }
        }
        adicionaElemento(s, bstPosA);
        adicionaElemento(s, bstPosB);
        s.numParCob = bstFO;
    }
    /*
    while ((s.numConIns > maxContReal_) || (s.numFaiCob > maxFaix_)) {
        bstFO = bstPosA = -1;

        for (int i = 0; i < s.numAreCom; i++) {
            foAntes = s.numParCob;
            removeElemento(s, i);
            calcParCob(s);

            if (s.numParCob > bstFO) {
                bstFO = s.numParCob;
                bstPosA = i;
            }

            adicionaElemento(s, i);
            s.numParCob = foAntes;
        }
        removeElemento(s, bstPosA);
        s.numParCob = bstFO;
    }
    */
}

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

void reconstruirSolucaoTri(Solucao& s) {
    double x;
    int foAntes = s.numParCob, bstFO, bstPosA, bstPosB, bstPosC;

    while ((s.numConIns < maxContReal_) && (s.numFaiCob < maxFaix_)) {
        bstFO = bstPosA = bstPosB = bstPosC = -1;
        for (int cont = 0; cont < s.numAreSem - 2; cont++) {
            for (int cont1 = cont + 1; cont1 < s.numAreSem - 1; cont1++) {
                for (int cont2 = cont1 + 1; cont2 < s.numAreSem; cont2++) {
                    adicionaElemento(s, cont);
                    adicionaElemento(s, cont1);
                    adicionaElemento(s, cont2);
                    calcParCob(s);

                    if (s.numParCob > bstFO) {
                        bstFO = s.numParCob;
                        bstPosA = cont;
                        bstPosB = cont1;
                        bstPosC = cont2;
                    }

                    removeElemento(s, s.numAreCom - 1);
                    removeElemento(s, s.numAreCom - 1);
                    removeElemento(s, s.numAreCom - 1);
                }
            }
        }
        adicionaElemento(s, bstPosA);
        adicionaElemento(s, bstPosB);
        adicionaElemento(s, bstPosC);
        s.numParCob = bstFO;
    }
    
    while ((s.numConIns > maxContReal_) || (s.numFaiCob > maxFaix_)) {
        removeElemento(s, rand() % s.numAreCom);
        calcParCob(s);
    }
}

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

void reconstruirSolucaoGul(Solucao& s) {
    int foAntes, bstFO, bstPos;

    while ((s.numConIns < maxContReal_) && (s.numFaiCob < maxFaix_)) {
        bstFO = bstPos = -1;
        for (int i = 0; i < s.numAreSem; i++) {
            adicionaElemento(s, i);
            calcParCob(s);

            if (s.numParCob > bstFO) {
                bstFO = s.numParCob;
                bstPos = i;
            }

            removeElemento(s, s.numAreCom - 1);
        }
        adicionaElemento(s, bstPos);
        s.numParCob = bstFO;
    }
}

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

void reconstruirSolucaoPior(Solucao& s, int aux) {
    adicionaElemento(s, s.numAreSem - aux);
    calcParCob(s);

    if (rand() % 2 == 0) 
        reconstruirSolucaoGul(s);
    else
        reconstruirSolucaoAle(s);
}

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

int selecionarHeuristicaDesconstrutiva() {
    int i;
    double sum, x, aux;
    sum = 0.0;
    for (i = 0; i < NUM_DES; i++)
        sum += pesoDesconstrutiva[i];
    x = rand() / (RAND_MAX + 1.0f);
    aux = 0.0;
    for (i = 0; i < NUM_DES; i++)
    {
        aux += (pesoDesconstrutiva[i] / sum);
        if (x < aux)
            return i;
    }
    return -1;
}

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

int selecionarHeuristicaReconstrutiva() {
    int i;
    double sum, x, aux;
    sum = 0.0;
    for (i = 0; i < NUM_REC; i++)
        sum += pesoReconstrutiva[i];
    x = rand() / (RAND_MAX + 1.0f);
    aux = 0.0;
    for (i = 0; i < NUM_REC; i++) {
        aux += (pesoReconstrutiva[i] / sum);
        if (x < aux) return i;
    }
    return -1;
}

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
/*
Solucao execALNS(Solucao& s) {
    Solucao sA, sB;
    clock_t hI, hF;

    int j = 0, heu_Des, heu_Rec, itMax, cnt;
    itMax = MAX_ITERACOES * numAre_;

    double x, T;

    memset(pesoDesconstrutiva, 1, sizeof(pesoDesconstrutiva));
    memset(pesoReconstrutiva, 1, sizeof(pesoReconstrutiva));
    memset(scoreDesconstrutiva, 0, sizeof(scoreDesconstrutiva));
    memset(scoreReconstrutiva, 0, sizeof(scoreReconstrutiva));
    memset(aparicoesDesconstrutiva, 0, sizeof(aparicoesDesconstrutiva));
    memset(aparicoesReconstrutiva, 0, sizeof(aparicoesReconstrutiva));

    hI = clock();
    criarSolucao(s);
    hF = clock();
    bstTime_ = ((double)hF - hI) / CLOCKS_PER_SEC;
    memcpy(&sA, &s, sizeof(s));
    foIni_ = s.numParCob;
    printf("\n>>> Sol. Ini.: %d\tTempo: %.2f\n\n", foIni_, bstTime_);

    solAva_ = 1;
    excTime_ = 0.0;

    while (excTime_ < MAX_TIME) {
        T = TEMP_INICIAL;

        while (T > TEMP_CONGELAMENTO) {
            for (int i = 0; i < itMax; i++) {
                solAva_++;
                i++; j++;
                memcpy(&sB, &s, sizeof(s));

                cnt = 2 + rand() % (int)floor(PER_REM*(s.numAreCom)) - 1;

                heu_Des = selecionarHeuristicaDesconstrutiva();
                heu_Rec = selecionarHeuristicaReconstrutiva();
                
                if (heu_Des == 0)
                    piorRemocao(sB, cnt);
                else if (heu_Des == 1)
                    shawRemocao(sB, cnt);
                else
                    aleRemocao(sB, cnt);
                aparicoesDesconstrutiva[heu_Des]++;
                
                if (heu_Rec == 0)
                    reconstruirSolucaoGul(sB);
                else if (heu_Rec == 1)
                    reconstruirSolucaoPior(sB);
                else
                    reconstruirSolucaoAle(sB);
                aparicoesReconstrutiva[heu_Rec]++;

                if (sB.numParCob > s.numParCob) {
                    memcpy(&s, &sB, sizeof(sB));

                    if (sB.numParCob > sA.numParCob) {
                        memcpy(&sA, &sB, sizeof(sB));

                        scoreDesconstrutiva[heu_Des] += S1;
                        scoreReconstrutiva[heu_Des] += S1;

                        hF = clock();
                        bstTime_ = ((double)hF - hI) / CLOCKS_PER_SEC;
                        printf(">>> Best Sol: %d\tTempo: %.2f\n", s.numParCob, bstTime_);
                    }
                    else {
                        scoreDesconstrutiva[heu_Des] += S2;
                        scoreReconstrutiva[heu_Des] += S2;
                    }
                }
                else {
                    x = (double)rand() / (double)RAND_MAX;

                    if (x < expl(-(sA.numParCob - sB.numParCob) / T)) {
                        memcpy(&s, &sB, sizeof(sB));
                        scoreDesconstrutiva[heu_Des] += S3;
                        scoreReconstrutiva[heu_Des] += S3;
                    }
                }

                if (j = SEG) {
                    j = 0;

                    for (int cont = 0; cont < NUM_DES; cont++) {
                        if (aparicoesDesconstrutiva[cont] != 0) {
                            pesoDesconstrutiva[cont] = (1 - FAT_REACAO) * pesoDesconstrutiva[cont] + (FAT_REACAO * scoreDesconstrutiva[cont]) / aparicoesDesconstrutiva[cont];
                            scoreDesconstrutiva[cont] = 0;
                        }
                    }

                    for (int cont = 0; cont < NUM_REC; cont++) {
                        if (aparicoesReconstrutiva[cont] != 0) {
                            pesoReconstrutiva[cont] = (1 - FAT_REACAO) * pesoReconstrutiva[cont] + (FAT_REACAO * scoreReconstrutiva[cont]) / aparicoesReconstrutiva[cont];
                            scoreReconstrutiva[cont] = 0;
                        }
                    }
                }
                
                hF = clock();
                excTime_ = ((double)hF - hI) / CLOCKS_PER_SEC;
                if (excTime_ >= MAX_TIME)
                    goto FIM;
            }
            T = TAXA_RESFRIAMENTO * T;
        }
    }

    hF = clock();
    FIM:;
    excTime_ = ((double)hF - hI) / CLOCKS_PER_SEC;
    foFin_ = sA.numParCob;
    printf("\n>>> Sol. Fin.: %d\tTempo: %.2f\tSol. Ava.: %d\n", foFin_, excTime_, solAva_);

    memcpy(&s, &sA, sizeof(sA));
    return sA;
}*/

//------------------------------------------------------------------------------

//==============================================================================

Solucao execALNS(Solucao& s) {
    Solucao sA, sB;
    clock_t hI, hF;

    int heu_Des, heu_Rec, itMax, cnt, iter, seg = 0;
    itMax = MAX_ITERACOES * numAre_;

    double x, T = 0;

    memset(pesoDesconstrutiva, 1, sizeof(pesoDesconstrutiva));
    memset(pesoReconstrutiva, 1, sizeof(pesoReconstrutiva));
    memset(scoreDesconstrutiva, 0, sizeof(scoreDesconstrutiva));
    memset(scoreReconstrutiva, 0, sizeof(scoreReconstrutiva));
    memset(aparicoesDesconstrutiva, 0, sizeof(aparicoesDesconstrutiva));
    memset(aparicoesReconstrutiva, 0, sizeof(aparicoesReconstrutiva));

    hI = clock();
    criarSolucao(s);
    hF = clock();
    bstTime_ = ((double)hF - hI) / CLOCKS_PER_SEC;
    memcpy(&sA, &s, sizeof(s));
    foIni_ = s.numParCob;
    printf("\n>>> Sol. Ini.: %d\tTempo: %.2f\n\n", foIni_, bstTime_);

    solAva_ = 1;
    excTime_ = 0.0;

    while (excTime_ < MAX_TIME) {
        if (T < TEMP_CONGELAMENTO)
            T = TEMP_INICIAL;
        iter = 0;

        while (iter < MAX_ITERACOES) {
            iter++;
            seg++;
            solAva_++;
            memcpy(&sB, &s, sizeof(s));

            cnt = 2 + rand() % (int)floor(PER_REM * (s.numAreCom)) - 1;

            heu_Des = selecionarHeuristicaDesconstrutiva();
            heu_Rec = selecionarHeuristicaReconstrutiva();

            if (heu_Des == 0)
                piorRemocao(sB, cnt);
            else if (heu_Des == 1)
                shawRemocao(sB, cnt);
            else
                aleRemocao(sB, cnt);
            aparicoesDesconstrutiva[heu_Des]++;

            if (heu_Rec == 0)
                reconstruirSolucaoGul(sB);
            else if (heu_Rec == 1)
                reconstruirSolucaoPior(sB, cnt);
            else
                reconstruirSolucaoAle(sB);
            aparicoesReconstrutiva[heu_Rec]++;

            if (sB.numParCob > s.numParCob) {
                memcpy(&s, &sB, sizeof(sB));

                if (sB.numParCob > sA.numParCob) {
                    memcpy(&sA, &sB, sizeof(sB));

                    scoreDesconstrutiva[heu_Des] += S1;
                    scoreReconstrutiva[heu_Des] += S1;

                    hF = clock();
                    bstTime_ = ((double)hF - hI) / CLOCKS_PER_SEC;
                    printf(">>> Best Sol: %d\tTempo: %.2f\n", s.numParCob, bstTime_);
                }
                else {
                    scoreDesconstrutiva[heu_Des] += S2;
                    scoreReconstrutiva[heu_Des] += S2;
                }
            }
            else {
                x = (double)rand() / (double)RAND_MAX;

                if (x < expl(-(sA.numParCob - sB.numParCob) / T)) {
                    memcpy(&s, &sB, sizeof(sB));
                    scoreDesconstrutiva[heu_Des] += S3;
                    scoreReconstrutiva[heu_Des] += S3;
                }
            }

            if (seg = SEG) {
                seg = 0;

                for (int cont = 0; cont < NUM_DES; cont++) {
                    if (aparicoesDesconstrutiva[cont] != 0) {
                        pesoDesconstrutiva[cont] = (1 - FAT_REACAO) * pesoDesconstrutiva[cont] + (FAT_REACAO * scoreDesconstrutiva[cont]) / aparicoesDesconstrutiva[cont];
                        scoreDesconstrutiva[cont] = 0;
                    }
                }

                for (int cont = 0; cont < NUM_REC; cont++) {
                    if (aparicoesReconstrutiva[cont] != 0) {
                        pesoReconstrutiva[cont] = (1 - FAT_REACAO) * pesoReconstrutiva[cont] + (FAT_REACAO * scoreReconstrutiva[cont]) / aparicoesReconstrutiva[cont];
                        scoreReconstrutiva[cont] = 0;
                    }
                }
            }

            hF = clock();
            excTime_ = ((double)hF - hI) / CLOCKS_PER_SEC;
            if (excTime_ >= MAX_TIME)
                goto FIM;
            T = TAXA_RESFRIAMENTO * T;
        }
    }

    hF = clock();
FIM:;
    excTime_ = ((double)hF - hI) / CLOCKS_PER_SEC;
    foFin_ = sA.numParCob;
    printf("\n>>> Sol. Fin.: %d\tTempo: %.2f\tSol. Ava.: %d\n", foFin_, excTime_, solAva_);
    printf("\nGul: %d  Pior: %d  Ale: %d\n", aparicoesReconstrutiva[0], aparicoesReconstrutiva[1], aparicoesReconstrutiva[2]);

    memcpy(&s, &sA, sizeof(sA));
    return sA;
}

//------------------------------------------------------------------------------

//==============================================================================