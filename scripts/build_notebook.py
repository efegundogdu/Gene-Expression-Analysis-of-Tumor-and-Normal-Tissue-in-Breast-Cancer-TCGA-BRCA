"""Build notebooks/analysis.ipynb (Turkish narrative, code calls src/ functions).

Run:  .venv/Scripts/python -I scripts/build_notebook.py
Then: .venv/Scripts/jupyter nbconvert --to notebook --execute --inplace notebooks/analysis.ipynb
"""

from pathlib import Path

import nbformat as nbf

ROOT = Path(__file__).resolve().parents[1]
nb = nbf.v4.new_notebook()
cells = []


def md(text):
    cells.append(nbf.v4.new_markdown_cell(text.strip()))


def code(text):
    cells.append(nbf.v4.new_code_cell(text.strip()))


md("""
# TCGA-BRCA: tümör ve normal doku arasında diferansiyel ekspresyon
## Batch etkisi ve yöntem seçiminin yanlış pozitiflere etkisi

Bu notebook projenin tüm aşamalarını `src/` içindeki fonksiyonlarla baştan sona çalıştırır.

- **Hızlı adımlar yeniden hesaplanır:** veri yükleme, ön işleme, rank testleri, benchmark ve zenginleştirme.
- **Ağır adımlar diskten yüklenir:** PyDESeq2, ComBat-seq, 2216 permütasyon çalışması ve iç içe çapraz doğrulama. Bu sonuçlar `src/de_methods.py`, `src/permutation.py` ve `src/downstream.py` ile üretildi; nasıl yeniden üretileceği README'de anlatılıyor.
- **Tekrarlanabilirlik kontrolü:** Yeniden hesaplanan sonuçlar, diskteki kayıtlı sonuçlarla karşılaştırılır.

Tüm rastgele işlemlerde seed = 42 kullanıldı.
""")

code("""
import sys, json, warnings, subprocess
from pathlib import Path
import numpy as np
import pandas as pd
from IPython.display import Markdown, display

warnings.filterwarnings("ignore")
ROOT = Path.cwd().resolve().parent
sys.path.insert(0, str(ROOT / "src"))
pd.set_option("display.width", 200)
pd.set_option("display.max_columns", 40)

import load, batch, de_methods as dm, permutation as pm, benchmark as bm, downstream as ds, plots, report

TABLES = ROOT / "results" / "tables"
RERUN_HEAVY = False  # True: rerun PyDESeq2 / ComBat-seq for Stage 3 (~25 min)
print(sys.executable)
""")

# ------------------------------------------------------------------ stage 1
md("""
## Aşama 1 — Veri yükleme ve etiketleme

GDC'den indirilen 1231 STAR "augmented gene counts" dosyasından `unstranded` sütunu okunur ve gen × örnek tamsayı matrisi oluşturulur. Baştaki `N_*` özet satırları atılır.

Aliquot barkodları `metadata.repository*.json` dosyasından alınır. Barkoddan şu bilgiler çıkarılır: TSS, hasta, örnek tipi, plate ve center.

Analizde sadece primer tümör (01) ve solid tissue normal (11) örnekleri kalır. Bir hastanın birden fazla tümör aliquot'u varsa alfabetik olarak ilk olanı tutulur.
""")

code("""
load.main()
samples = pd.read_csv(ROOT / "data/processed/samples.csv")
samples["condition"].value_counts()
""")

md("""
**Yorum.** Sonuçlar beklenen sayılarla örtüşüyor:
- Dosyaların dağılımı: 1111 primer tümör, 113 normal ve 7 metastaz (metastazlar çıkarıldı).
- 11 hastanın birden fazla tümör aliquot'u vardı. Tekrarlar atılınca 1095 tümör ve 113 normal kaldı.
- 113 hastanın hem tümör hem normal örneği var. Bunların 89'unda iki örnek aynı plate'te.
- Plate sayısı tekrarlar atılmadan önce 45, sonra 44. A277 plate'inde sadece tekrar aliquot'lar vardı.
- Center tüm örneklerde 07 olduğu için batch değişkeni olarak kullanılmıyor.
""")

# ------------------------------------------------------------------ stage 2
md("""
## Aşama 2 — Ön işleme ve batch tanılaması

- **Gen filtresi:** Protein-coding genler alınır. Bir gen, en az 113 örnekte (normal örnek sayısı) CPM ≥ 1 ise tutulur.
- **Normalizasyon:** DESeq2 median-of-ratios size factor'leri kullanılır.
- **Aykırı örnekler:** PCA ve kütüphane boyutuna göre işaretlenir ama silinmez.
- **Batch etkisi:** İlk 10 PC'nin ne kadarını plate, TSS ve condition'ın açıkladığı R² ile ölçülür.
- **Plate gruplama:** 10'dan az örnekli plate'ler, işlem sırasına göre komşu plate'lerle birleştirilir; sonuç `batch_group` sütunudur.
""")

code("""
res2 = batch.main()
plots.make_stage2_figures(res2)
print("genes kept:", res2["filt"].shape[0], "| outliers flagged:", int(res2["flags"]["outlier_flag"].sum()))
t = res2["test_plate"]
print(f"plate x condition: chi2={t['chi2']:.1f}, p={t['p']:.2g}, Cramer's V={t['cramers_v']:.3f}")
print("batch groups:", len(res2["batch_table"]), "| with both conditions:", len(res2["both"]),
      "| samples in them:", int(res2["samples"]["batch_has_both"].sum()))
res2["pcv"].round(3)
""")

md("""
![Fig 2](../results/figures/fig2_pca.png)

![Fig 3](../results/figures/fig3_plate_condition_heatmap.png)

**Yorum.**
- **PC1:** Tümör ile normali ayırıyor (condition R² = 0.50).
- **PC3:** Neredeyse yarısını plate açıklıyor (R² = 0.50). Bu etki işlem sırasına bağlı: geç plate'ler sadece tümör içeriyor ve PC3'te pozitif tarafta (`figS_pc3_by_plate.png`).
- **Plate ile condition güçlü şekilde iç içe:** Cramér's V = 0.55.
  - A089 plate'inde 22 normal var, hiç tümör yok.
  - 27 plate'te sadece tümör var.
  - 38 batch grubunun yalnızca 15'i iki koşulu birden içeriyor. Bu 15 grupta 615 örnek var.
""")

# ------------------------------------------------------------------ stage 3
md("""
## Aşama 3 — Yedi DE protokolü

| | Yöntem |
|---|---|
| P1 | PyDESeq2 `~condition` (apeGLM ile shrink edilmiş log2FC) |
| P2 | PyDESeq2 `~batch_group + condition` |
| P3 | ComBat-seq (R sva), ardından PyDESeq2 `~condition` |
| P4 | Normalize count'larda Wilcoxon rank-sum |
| P5 | ComBat-seq ile düzeltilmiş count'larda Wilcoxon |
| P6 | van Elteren stratifiye Wilcoxon (strata = `batch_group`) |
| P7 | Koşul içi 95. persentilde winsorize edilmiş count'larla PyDESeq2 |

DEG tanımı: FDR (Benjamini-Hochberg) < 0.05 ve |log2FC| > 1.

Önce van Elteren testinin doğrulama testleri çalıştırılır. Ardından hızlı protokoller (P4, P6) yeniden hesaplanıp kayıtlı tablolarla karşılaştırılır.
""")

code("""
out = subprocess.run([sys.executable, "-m", "pytest", "-q", str(ROOT / "tests")], capture_output=True, text=True)
print(out.stdout[-300:])
""")

code("""
counts, samples_i, genes = dm.load_data()
check = {}
for p in ["P4", "P6"]:
    new, _ = dm.run_protocol(p, counts, samples_i)
    old = dm.load_protocol(p)
    check[p] = {"same DEG set": set(new.index[new.is_DEG]) == set(old.index[old.is_DEG]),
                "max |p diff|": float(np.nanmax(np.abs(new["pvalue"].to_numpy() - old["pvalue"].to_numpy())))}
if RERUN_HEAVY:
    dm.main(["P1", "P2", "P3", "P5", "P7", "--force"])
pd.DataFrame(check).T
""")

code("""
r3 = dm.summarize_protocols()
plots.fig4_deg_counts(r3["info"]); plots.fig_deg_jaccard(r3["jaccard"])
display(r3["info"][["n_used", "n_DEG", "n_up", "n_down", "seconds"]].round(1))
common = dm.common_lfc_sensitivity()
display(Markdown("**DEG Jaccard: own log2FC vs common log2FC**"))
display(pd.concat({"own": r3["jaccard"].round(2), "common": common["jaccard"].round(2)}, axis=1))
""")

md("""
![Fig 4](../results/figures/fig4_deg_counts.png)

**Yorum.**
- **Markerlar:** Yedi protokolün hepsi beş kanonik markerı doğru yönde buluyor.
- **DEG sayısını log2FC eşiği belirliyor:** 1208 örnekle genlerin %80'inden fazlası FDR < 0.05 çıkıyor.
- **Ortak log2FC'nin etkisi:** Kendi log2FC'leriyle DESeq2 ve rank aileleri belirgin şekilde ayrışıyor (Jaccard yaklaşık 0.6). Ortak bir DESeq2 log2FC kullanılınca P1 ile P4 arasındaki Jaccard 0.59'dan 0.94'e çıkıyor.
- **Sonuç:** Ailelerin ayrışmasının çoğu testten değil, log2FC'nin ortalama ya da medyanla tanımlanmasından geliyor.
- **P7 notu:** P7 başta tüm örnekler üzerinden winsorize edilmişti. Bu, %9'luk normal grubu kırptığı için ADH1B'nin log2FC'sini −3.6'ya düşürdü; o yüzden winsorization koşul içinde yapıldı.
""")

# ------------------------------------------------------------------ stage 4
md("""
## Aşama 4 — Permütasyon analizi (projenin özgün katkısı)

Bu bölüm 2216 protokol çalışmasının diskteki sonuçlarını yükler.

- **Null A (yöntem kaynaklı):** Etiketler tüm örnekler arasında karıştırılır.
- **Null B (batch kaynaklı):** Sadece tümör örnekleri kullanılır, dolayısıyla bulunan her DEG yanlış pozitiftir. Sahte normaller iki şekilde seçilir:
  - **B-rand:** Rastgele seçilir.
  - **B-conf:** Gerçek verideki plate yapısını taklit eder: yaklaşık %20'si tamamen sahte-normal plate'lerde, plate × etiket Cramér's V değeri gerçek veriye yakın.
  - **FP_batch** = FP_conf − FP_rand. İkisi aynı replikasyonda, aynı örnekler ve aynı sayıda sahte normalle hesaplanır.
- **Null C (ikincil):** Etiketler sadece iki koşulu da içeren 15 batch grubu içinde karıştırılır.
- **Hangi veri:** Rank testleri (P4, P6) tam veride B = 100 ile çalışır. Diğer protokoller alt örneklemde B = 50 ile çalışır: Null A'da 113 normal + 339 tümör, Null B'de 452 tümör.
""")

code("""
r4 = pm.run_summary()
cols = ["protocol", "FP_A", "FP_rand", "FP_conf", "FP_batch", "FP_batch_ci_lo", "FP_batch_ci_hi", "FP_C", "FDR_A"]
for thr, d in r4["summary"].groupby("threshold", sort=False):
    display(Markdown(f"**{thr}**")); display(d[cols].round(2).set_index("protocol"))
display(pd.read_csv(TABLES / "perm_label_diagnostics.csv", index_col=0).round(3))
""")

md("""
![Fig 5](../results/figures/fig5_false_positives.png)

**Yorum.**
- **Li et al. 2022 tekrarlandı:** DESeq2 tabanlı protokoller, etiketler rastgele karıştırıldığında bile replikasyon başına yüzlerce yanlış pozitif veriyor (P1: 357). Rank testleri neredeyse sıfır.
- **Ama tahmini FDR düşük:** Gerçek veride yaklaşık 12 bin gen anlamlı olduğu için tahmini FDR %3'ün altında kalıyor.
- **Batch karışıklığı P4 ve P5'i bozuyor:** Plate'lerle karışık etiketlerde düz Wilcoxon (P4) ortalama 378, ComBat-seq + Wilcoxon (P5) 374 batch kaynaklı yanlış pozitif veriyor.
- **ComBat-seq korumuyor, zarar veriyor:** Etiketler batch'le karışık olduğunda ComBat-seq, koşul kovaryatını korurken batch farkını koşul sinyaline aktarıyor (Nygaard et al. 2016). P5'in bu yanlış pozitiflerinin çoğu |log2FC| > 1 eşiğini de geçiyor (240).
- **van Elteren (P6) her durumda dayanıklı:** FP_batch = 0.2.
- **Null C'ye dikkat:** Etiketler gerçek etiketlerle hâlâ hafifçe ilişkili (korelasyon 0.07–0.11). Yanlış pozitif sayısı bu korelasyonla artıyor; Spearman P1'de 0.83. Yani Null C batch ile biyolojiyi karıştırıyor.
""")

code("""
bs = pd.read_csv(TABLES / "batch_sensitive_genes.csv", index_col=0)
pd.DataFrame({p: {"batch-sensitive": int(bs[f"batch_sensitive_{p}"].sum()),
                  "batch-specific": int(bs[f"batch_specific_{p}"].sum())} for p in pm.PROTOCOL_ORDER})
""")

md("""
**P7 duyarlılık analizi (global winsorization, "P7g").** Yang et al. 2024'ün prosedürü uygulandı, ama persentil koşul içinde değil tüm örnekler üzerinden hesaplandı. Null A alt örnekleminde 50 replikasyonla çalıştırıldı.
""")

code("""
pd.read_csv(TABLES / "p7_winsorization_comparison.csv", index_col=0).round(3)
""")

md("""
**Yorum.**
- **P7g:** Null A'da yanlış pozitifleri P1'e göre %94 azaltıyor. |log2FC| > 1 eşiği eklenince azalma %98.7. Bu, Yang et al.'ın bildirdiği %98.2 ile uyumlu.
- **Koşul içi P7:** Azalma sadece %3.
- **Bedel:** Global winsorization gerçek veride daha az DEG buluyor (3401'e karşı 3785).
- **Ana tablolar:** P7 ana tablolarda koşul içi versiyon olarak kaldı; P7g burada duyarlılık analizi olarak raporlanıyor.
""")

# ------------------------------------------------------------------ stage 5
md("""
## Aşama 5 — Benchmark

**Referans:** Aynı plate'te eşleşmiş 89 hasta (178 örnek) kullanıldı. Bu çiftlerde tümör ile normal aynı plate'te olduğu için plate karışıklığı yok. İki yöntemle DEG bulundu:
- Hasta içi farklar üzerinde Wilcoxon signed-rank.
- PyDESeq2 `~patient + condition`.

Referans listesi iki yöntemin kesişimi: ikisinde de aynı yönde DEG olan genler.
""")

code("""
r5 = bm.run_all()
print(json.dumps(r5["reference_info"], indent=1))
display(pd.concat({"own log2FC": r5["benchmark"][["precision", "recall", "jaccard"]],
                   "common log2FC": r5["benchmark_common_lfc"][["precision", "recall", "jaccard"]]}, axis=1).round(3))
""")

md("""
![Fig 6](../results/figures/fig6_upset.png)

![Fig 7](../results/figures/fig7_markers.png)

**Yorum.**
- **Referans:** 3262 gen (1137 up, 2125 down). İki yöntem arasında yön çelişkisi yok.
- **Kendi log2FC'leriyle:** Rank protokolleri referansa en yakın (Jaccard yaklaşık 0.81).
- **Ortak log2FC ile:** Bu üstünlük kayboluyor; P4 0.59, P1 0.57. Öne çıkanlar batch'i ele alan protokoller: P2, P3, P5 ve P6, Jaccard 0.66–0.69.
- **Neden:** Referansın Wilcoxon yarısı da medyana dayalı log2FC kullanıyor, bu da rank protokollerini kayırıyor.
- **Bağımsızlık sınırı:** Referanstaki 178 örnek protokollerin kullandığı tam veride de var; referans bağımsız bir veri seti değil.
- **Kesişimler (Fig 6):** 7 protokolün ortak kesişimi 2452 gen. Sonraki en büyük kesişimler aileye özgü: 491 gen sadece DESeq2 ailesinde, 395 gen sadece rank ailesinde.
""")

# ------------------------------------------------------------------ stage 6
md("""
## Aşama 6 — Downstream etki

- **Zenginleştirme:** GO BP 2023 ve KEGG 2021 kütüphaneleri kullanıldı. Arka plan, test edilen 15317 gen. Her protokolün ilk 20 teriminin protokoller arası Jaccard benzerliği hesaplandı.
- **Sınıflandırma:** İç içe (nested) 5 katlı çapraz doğrulama. Gen seçimi her eğitim katmanında yeniden yapılıyor, böylece test katmanından bilgi sızmıyor. Sonuçlar diskten yüklenir.
- **Robust core imzası:** P2, P3, P5 ve P6'nın dördünde de aynı yönde DEG olan ve B-conf replikasyonlarının %5'inden azında anlamlı çıkan genler.
""")

code("""
r6 = ds.run_all()
display(r6["enrichment"]["summary"].pivot_table(index="protocol", columns=["library", "direction"], values="n_terms_FDR<0.05"))
display(r6["classification"][["name", "model", "n_genes", "roc_auc", "f1_normal", "stability_selected", "stability_L1-LR", "stability_RF"]].round(4))
for n in ["robust_core", "robust_core_v2"]:
    print(n, json.loads((TABLES / f"{n}_info.json").read_text()))
""")

md("""
![Fig 8](../results/figures/fig8_term_jaccard.png)

![Fig 9](../results/figures/fig9_classification.png)

**Yorum.**
- **Zenginleştirme:** Biyoloji tutarlı: up genlerde mitoz ve hücre döngüsü, down genlerde anjiyogenez, kas ve sinir sistemi terimleri. Ama ilk 20 terimin örtüşmesi yöntem ailesine göre değişiyor; GO BP down'da DESeq2 ile rank aileleri arasında Jaccard 0.29–0.43.
- **Sınıflandırma:** Tümör ile normal bu veride çok kolay ayrıldığı için AUC tüm protokollerde ≥ 0.9994. Protokol seçimi sınıflandırma başarısını değiştirmiyor.
- **Gen kararlılığı:** Seçilen DEG listeleri katmanlar arasında kararlı (Jaccard yaklaşık 0.92). L1 modelinin gerçekten kullandığı genler ise değişken (0.18–0.31).
- **Robust core v1:** 1664 gen. Elenen genlerin çoğunu P5'in kendi batch kaynaklı yanlış pozitifleri belirliyor; ADH1B, FABP4 ve LEP bu yüzden düşüyor.
- **Robust core v2:** Sıklık filtresi sadece P2, P3 ve P6'ya uygulanıyor. Sonuç 2010 gen; ADH1B ve FABP4 geri dönüyor.
- **İmza karşılaştırması:** Her iki imza, P1'in tam DEG listesinin yaklaşık yarısı kadar genle aynı başarıyı veriyor. Bu karşılaştırma nested değil, çünkü listeler tüm örneklerle seçildi.
""")

# ------------------------------------------------------------------ stage 7
md("""
## Aşama 7 — Özet tablo (Tablo 1)
""")

code("""
t1 = report.build_table1()
t1.drop(columns=["method"]).round(3).T
""")

md("""
## Sonuç

1. **Tasarım batch'le iç içe:** TCGA-BRCA'da plate ile tümör/normal durumu güçlü şekilde iç içe (Cramér's V = 0.55). Örneklerin yarısı tek koşullu batch'lerde.
2. **Rank testleri ve batch:** Rastgele etiketlerde rank testleri yanlış pozitif vermiyor. Ama etiketler plate'lerle karışık olduğunda düz Wilcoxon yüzlerce yanlış pozitif veriyor. Stratifiye van Elteren testi bu karışıklığa karşı tamamen dayanıklı.
3. **ComBat-seq bu tasarımda zararlı:** Koşul kovaryatıyla ComBat-seq uygulamak batch'i düzeltmiyor, aksine batch farkını koşul sinyaline aktarıyor.
4. **DESeq2 ve winsorization:** DESeq2'nin yöntem kaynaklı yanlış pozitifleri gerçek, ama bu veride tahmini FDR %3'ün altında. Global winsorization bu yanlış pozitiflerin %94'ünü gideriyor.
5. **Protokol farkının kaynağı:** Protokoller arasındaki DEG farkının büyük kısmı log2FC tanımından geliyor. Sınıflandırma başarısı protokolden etkilenmiyor; zenginleştirmede ise ilk terimler yöntem ailesine göre değişiyor.
""")

nb["cells"] = cells
nb["metadata"]["kernelspec"] = {"name": "python3", "display_name": "Python 3", "language": "python"}
nb["metadata"]["language_info"] = {"name": "python"}
out = ROOT / "notebooks" / "analysis.ipynb"
out.parent.mkdir(parents=True, exist_ok=True)
nbf.write(nb, out)
print("wrote", out)
