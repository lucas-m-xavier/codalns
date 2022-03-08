#ifndef uPCTFH
#define uPCTFH

#include <iostream>

#define MAX_NOS 2500 // n�mero m�ximo de n�s - baseado no tamanho das inst�ncias
#define MAX_ARE 2600 // n�mero m�ximo de arestas - baseado no tamanho das inst�ncias
#define MAX_PAR 900  // n�mero m�ximo de pares O-D - baseado no tamanho das inst�ncias
#define NUM_DES 3 // n�mero de heuristicas desconstrutivas
#define NUM_REC 3 // n�mero de heuristicas reconstrutivas

// ------------------------------- ESTRUTURAS -------------------------------
typedef struct tAresta
{
 int id;   // identificador real da aresta
 int no1;  // id (posi��o no vetor de n�s) do primeiro n� da aresta
 int no2;  // id (posi��o no vetor de n�s) do segundo n� da aresta
 int nCon; // n�mero de contadores da aresta
 int nFai; // n�mero de faixas da aresta
}Aresta;

typedef struct tSolucao
{
	//int vetAre[MAX_ARE]; // vetor bin�rio de arestas (1 - contador instalado; 0 - caso contr�rio)
	int vetIDCon[MAX_ARE]; // vetor de ID
	int vetIDSem[MAX_ARE]; // vetor de ID
	int numAreCom; // n�mero de arestas COM contador
	int numAreSem; // n�mero de arestas SEM contador
	int numConIns; // n�mero de contadores instalados
	int numFaiCob; // n�mero de faixas cobertas 
	int numParCob; // n�mero de pares cobertos
}Solucao;
//---------------------------------------------------------------------------


//---------------------------- VARI�VEIS GLOBAIS ----------------------------
// ----- Dados de entrada
int numPar_; // n�mero de pares O-D 
int numNos_; // n�mero de n�s na rede
int numAre_; // n�mero de arestas na rede
int numAreVir_; // n�mero de arestas virtuais
int numTotAre_; // n�mero total de arestas (inst�ncia completa)
int numTotFai_; // n�mero total de faixas das arestas (inst�ncia completa)
int *vetIdPar_; // vetor com o id real do n�s usados como pares O-D
int *vetIdNos_; // vetor com o id real dos n�s
Aresta *vetArestas_; // vetor de arestas
// ----- Rede
int *vetQtdAdj_;  // vetor com a quantidade de n�s (e arestas) adjacentes a cada n�
int **matNosAdj_; // matriz com os n�s adjacentes a cada n�
int **matAreAdj_; // matriz com as arestas adjacentes a cada n�
// ----- Resultados
int foIni_,foFin_;        // fun��o objetivo das solu��es inicial e final
double bstTime_,excTime_; // tempo para encontrar a melhor solu��o e tempo total
int solAva_;              // n�mero de solu��es avaliadas pelo CS
// ----- Vari�veis auxiliares
int maxPar_;      // n�mero m�ximo de pares cobertos
int maxCont_;     // n�mero m�ximo de contadores a serem instalados  
int maxFaix_;     // n�mero m�ximo de faixas a serem cobertas
int maxContReal_; // n�mero m�ximo REAL de contadores (definido com base no limite de contadores e faixas)
//-----ALNS
int cost[MAX_ARE];
int pos[MAX_ARE];
int pesoDesconstrutiva[NUM_DES];
int pesoReconstrutiva[NUM_REC];
int scoreDesconstrutiva[NUM_DES];
int scoreReconstrutiva[NUM_REC];
int aparicoesDesconstrutiva[NUM_DES];
int aparicoesReconstrutiva[NUM_REC];
static double cstRgtVec__[MAX_ARE];     // vetor com os custos utilizados na heurística de inserção k-regret (rgtInsHeu)
static int cntRgtVec__[MAX_ARE];        // vetor com os contadores utilizados na heurística de inserção k-regret (rgtInsHeu)
static int posRgtVec__[MAX_ARE];        // vetor com as posições utilizados na heurística de inserção k-regret (rgtInsHeu)

//---------------------------------------------------------------------------


//--------------------------------- M�TODOS ---------------------------------
// ------------ SA
void execSA(Solucao &s);
void gerVizinho(Solucao &s);
// ------------ Solu��o
void criarSolucao(Solucao &s);
void calcParCob(Solucao &s);
void calcParCob2(Solucao& s, const int* vetAre);
int calcPar(Solucao &s,const int &n1, const int *vetAre);
// ------------ Entrada e Sa�da
void lerInstancia(char *arq);
void escreverSolucao(Solucao &s,FILE *f);
void escreverResultado(Solucao &s,char *path);
// ------------ Auxiliares
void montarRede();
void limparMemoria();
//------------------ ALNS
void calculaCusto(Solucao& s);
void ordenaCusto(Solucao& s);
void ordenaCustoShaw(int aux, Solucao& s);
void adicionaElemento(Solucao& s, int aux);
void removeElemento(Solucao& s, int aux);
void piorRemocao(Solucao& s, int aux);
void aleRemocao(Solucao& s, int aux);
void shawRemocao(Solucao& s, int aux);
void reconstruirSolucaoAle(Solucao& s);
void reconstruirSolucaoPar(Solucao& s);
void reconstruirSolucaoTri(Solucao& s);
void reconstruirSolucaoGul(Solucao& s);
void reconstruirSolucaoPior(Solucao& s, int aux);
int selecionarHeuristicaDesconstrutiva();
int selecionarHeuristicaReconstrutiva();
Solucao execALNS(Solucao& s);
#endif