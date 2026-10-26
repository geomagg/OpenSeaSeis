#!/bin/bash
# Sintetico da analise de velocidade por focalizacao (SPK_XVSPMIG1 + SPK_TAU2VFIX)
# Uso: ./run_all.sh [NCPU]      (NCPU = modelagens em paralelo, padrao 2)
set -e
cd "$(dirname "$0")"
NCPU=${1:-2}
command -v seaseis >/dev/null || { echo "seaseis nao esta no PATH"; exit 1; }

echo "== 1. modelos de velocidade"
python3 mkmod.py 512 512 10 1.0 VTRUE10      # modelo verdadeiro para a modelagem (10 m)
python3 mkmod.py 512 200 25 1.0 VTRUE25      # verdadeiro, para a migracao (25 m)
python3 mkmod.py 512 200 25 0.9 V90_25       # 90% nos sedimentos (modelo inicial errado)

echo "== 2. modelagem dos 17 nos (~3-4 min por no)"
if [ ! -f node17.cseis ]; then seq 1 17 | xargs -P $NCPU -I{} ./runfd.sh {}; else echo "   node*.cseis ja existem, pulando"; fi
seaseis -D I=node9.cseis -D O=node9.sgy -f tosgy.flow -o tosgy.log > /dev/null

echo "== 3. migracao com a velocidade verdadeira (controle) e com 90%"
seaseis -D V=VTRUE25 -D P=TRUE -f MIG.flow -o mig_TRUE.log > /dev/null
seaseis -D V=V90_25  -D P=V90  -f MIG.flow -o mig_V90.log  > /dev/null
for p in TRUE V90; do seaseis -D P=$p -f T2V.flow -o t2v_$p.log > /dev/null; done

echo "== 4. picks -> Vint -> modelo novo, e iteracoes"
python3 update.py V90_vfix.sgy 9 VIT1_25
for k in 1 2 3; do
  seaseis -D V=VIT${k}_25 -D P=IT$k -f MIG.flow -o mig_IT$k.log > /dev/null
  seaseis -D P=IT$k -f T2V.flow -o t2v_IT$k.log > /dev/null
  echo "   -- picks da migracao $k:"
  python3 update.py IT${k}_vfix.sgy 9 VIT$((k+1))_25
done

echo "== 5. desvio tau e profundidades por migracao, figuras"
python3 stats.py
python3 figs.py
echo "pronto: fig_true.png fig_v90.png fig_it2.png fig_dados.png"
