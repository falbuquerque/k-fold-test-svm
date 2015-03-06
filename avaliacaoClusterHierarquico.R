###############################################################################
# Author: Felipe Albuquerque
###############################################################################

#
# Criação da classe que representa um k-grupo. Guarda uma matriz de dados no atributo "grupo"
#
setClass("k.grupo", representation(grupo = "matrix"))

#
# Valida o valor de k, que deve estar entre 1 e N. 
# Retorna 1 quando k é inválido e 0 quando é válido
#
isKValido = function(k, n) {((k > 1) & (k < n))}

#
# Divide os dados em k grupos
# 
# dados: os dados que devem ser divididos
# k: o número de grupos. Deve estar entre 1 e o número de dados
#
dividirEmKGrupos = function(dados, k) {
    
    # Validação de k
    linhas = length(dados[,1])
    
    if (!isKValido(k, linhas)) {
        stop(paste("O valor de k deve estar entre 1 e", linhas))
    }
    
    # Inicialização do retorno. Criação das estruturas dos grupos
    cols = length(dados[1,])
    conjuntoRetornado = c()
    
    for (j in 1:k) {
        conjuntoRetornado = c(conjuntoRetornado, new("k.grupo", grupo = matrix(nrow = 0, ncol = cols)))
    }
    
    # Criação dos grupos    
    
    for (i in 1:linhas) {
        pos = i %% k # Seleção do grupo baseada no resto da divisão
        
        if (pos == 0) { # Como o índice começa em "1", utilizar "k" ao invés de "0"
            pos = k
        }
        
        # Adiciona o dado ao grupo selecionado
        conjuntoRetornado[[pos]]@grupo = rbind(conjuntoRetornado[[pos]]@grupo, dados[i,])
    }
    
    conjuntoRetornado
}

#
# Retorna a matriz sem a coluna de classe
# 
# m: a matriz
# indiceColunaClasse: o índice da coluna de classe
#
extrairMatrizSemColunaDeClasse = function(m, indiceColunaClasse) {
    linhas = length(m[,1])
    cols = length(m[1,])
    
    if (indiceColunaClasse == cols) { # coluna de classe é a última
        m[1:linhas, (1:cols - 1)]
    } else {
        x = m[1:linhas, 1:(indiceColunaClasse -1)]
        cbind(m[1:linhas, (indiceColunaClasse + 1):cols])
    }
    
}

#
# Retorna a coluna de classe da matriz
# 
# m: a matriz
# indiceColunaClasse: o índice da coluna de classe
#
extrairColunaDeClasse = function(m, indiceColunaClasse) {
    linhas = length(m[,1])
    as.matrix(m[1:linhas, indiceColunaClasse:indiceColunaClasse], nrow = linhas, ncol = 1, byrow = TRUE)
}

#
# Avalia o classificador baseado em agrupamento hierárquico utilizando a técnica de cross-validation
#
avaliarAgrupamentoHierarquico = function(dados, k, indiceColunaClasse, funcaoMedicaoDistancia) {
    m = as.matrix(dados)
    linhas = length(m[,1])
    cols = length(m[1,])
    
    # Validação de k
    if (!isKValido(k, linhas)) {
        stop(paste("O valor de k deve estar entre 1 e", linhas))
    }
    
    # Execução da avaliação
    grupos = dividirEmKGrupos(m, k)
    
    print(paste("Validação cruzada para k =", k))
    print(paste("Função de medição de distância selecionada:", funcaoMedicaoDistancia))
    
    tp = 0
    fp = 0
    fn = 0
    tn = 0
    total = 0
    
    for (i in 1:k) {
        
        # Escolhe o grupo de testes
        grupoTeste = i
        conjuntoTestes = grupos[[grupoTeste]]@grupo
        
        # Como é não supervisionado, não há modelo de treinamento
    
        # Teste
        x = extrairMatrizSemColunaDeClasse(conjuntoTestes, indiceColunaClasse)
        y = extrairColunaDeClasse(conjuntoTestes, indiceColunaClasse)
        
        cluster = hclust(dist(x), funcaoMedicaoDistancia)
        
        # Algoritmo preparado apenas para classificação binária
        resultadoTestes = table(cutree(cluster, k = 2), y) 
        
        # Incremento dos valores de medidas
        tp = tp + resultadoTestes[1, 1]
        fp = fp + resultadoTestes[1, 2]
        fn = fn + resultadoTestes[2, 1]
        tn = tn + resultadoTestes[2, 2]
        total = total + length(conjuntoTestes[,1])
    }
    
    # Cálculos das medidas
    p = tp + fn
    n = tn + fp
    taxaFalsaAceitacao = fp / total
    taxaFalsaRejeicao = fn / total
    sensibilidade = tp / p
    especificidade = tn / n
    erroMedio = (fp + fn) / total
    
    print(paste("Positivos:", p))
    print(paste("Negativos:", n))
    print(paste("Verdadeiros Positivos:", tp))
    print(paste("Falsos Positivos:", fp))
    print(paste("Falsos Negativos:", fn))
    print(paste("Verdadeiros Negativos:", tn))
    print(paste("Erro médio:", erroMedio))
    print(paste("Taxa de Falsa Aceitação:", taxaFalsaAceitacao))
    print(paste("Taxa de Falsa Rejeição:", taxaFalsaRejeicao))
    print(paste("Sensibilidade:", sensibilidade))
    print(paste("Especificidade:", especificidade))
}

# Carga dos dados
f = file("D:\\USP\\PPgSI - Mestrado\\2010-2_Sem\\SIN 5007-1 - Reconhecimento de Padrões\\Atividades\\Atividade 6\\census_income_colunas-5_6-selecionadas.csv")
censusIncome = read.table(f, sep = ",", colClasses = c("integer", "integer", "integer"), header = TRUE)

require("graphics")

# Testes para diferentes funções de medição de distância 
funcoes = c("ward", "single", "complete", "average", "mcquitty", "centroid")

for (func in funcoes) {
    avaliarAgrupamentoHierarquico(dados = censusIncome, k = 10, indiceColunaClasse = 3, funcaoMedicaoDistancia = func)
}

plot(wardCluster)
plot(singleCluster)
plot(completeCluster)
plot(averageCluster)
plot(mcquittyCluster)
plot(centroidCluster)