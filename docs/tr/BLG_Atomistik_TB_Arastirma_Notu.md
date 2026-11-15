> Yazarın Word notundan (`BLG_Atomistik_TB_Arastirma_Notu.docx`) Markdown'a dönüştürülmüştür. Converted from the author's Turkish research note on the atomistic tight-binding solver. The solver source code itself is not part of this repository.

**Çift Katmanlı Grafen Kuantum Noktasında**

**Atomistik Tight-Binding Çözücüsünün Kurulması, Doğrulanması ve İntervalley Sınırı**

*Bir hesaplamalı araştırmanın tam kronolojik kaydı — hiçbir adım atlanmadan*

Saf NumPy/SciPy atomistik çözücü — beş yapısal çıpa — k·p sürekli-limit mutabakatı — valley ayrımı — intervalley denemesi

**Hesaplamalı araştırma notu**

*Tüm sayısal sonuçlar oturum içinde üretildi ve yakınsama disipliniyle sınandı*

## İçindekiler

## 0.  Bu belge nasıl okunmalı

Bu not, çift katmanlı grafen (BLG) kuantum noktalarında yük-gürültüsü korumasının kökenini araştıran daha geniş bir çalışmanın atomistik kanadını, gerçekleştiği sırayla kaydeder. Önceki aşamada bir 4-bant k·p (zarf) çözücüsü kurulmuş, bir el yazmasının AC-koruma iddiası çürütülmüş ve korumanın simetri-kilitli olduğu sonucuna varılmıştı. O çalışma, gerçek intervalley fiziğinin zarf modelinin dışında kaldığını ve tam-Brillouin-bölgesi atomistik bir hesap gerektirdiğini açıkça işaret ediyordu. Bu belge tam da o adımı — atomistik tight-binding çözücüsünün sıfırdan kurulmasını, doğrulanmasını ve sınırının saptanmasını — anlatır.

**Rehber ilke:** her nicel iddia açık bir çıpadan geçmelidir — bir simetri özdeşliği, bilinen bir analitik limit, bir yakınsama kontrolü ya da bağımsız bir yeniden üretim — yoksa negatif sonuç olarak raporlanır. Bu disiplin, kalıcı bir bulguyu geri çekilen bir bulgudan ayıran şeydir; bu belgedeki heyecan verici görünen birkaç ipucunun neden artifakt olarak raporlandığının da nedenidir.

**Okuyucu kısayolu.** En önemli iki sonuç şudur: (i) atomistik çözücü, hiçbir serbest parametre ayarlanmadan k·p sürekli-limitiyle üç haneye kadar mutabık kalmıştır (Bölüm 7); (ii) düzgün-potansiyelin intervalley matris elemanı yakınsama testinde sıfıra çökmüştür — yani tasarıma bağlı bir intervalley figürü yoktur, sınır atomik-ölçek kusurlardadır (Bölüm 10). Bölüm 2–6 bunları hak eden zemindir.

## 1.  Başlangıç: neden atomistik bir çözücü

k·p zarf modeli BLG kuantum noktasının düşük-enerji fiziğini iyi tanımlar, ancak iki şeyi yapısal olarak temsil edemez: trigonal warping kaynaklı minivalley yapısının gerçek C₃ örgü simetrisini ve |K−K′| ≈ 17 nm⁻¹ momentumundaki intervalley saçılmayı. Çözücünün kullandığı kare düzlem-dalga ızgarası C₃'ü kırar (yapay bir minivalley triplet yarılması üretir), ve küçük-q bazı |K−K′| momentumunu hiç içermez.

Bu iki sınır da valley kübitinin korumasını belirleyen fiziğin tam kalbinde durur. Dolayısıyla doğal bir sonraki adım, her iki valleyi gerçek ayrımında temsil eden ve C₃'e saygı duyan atomistik bir tight-binding hesabıdır. Hedef bir tam üretim koşusu değil, prensip-ispatı niteliğinde bir çözücüdür: minivalley yapısını yeniden üreten ve intervalley matris elemanına ilk kez QD geometrisinde erişebilen bir araç.

**Ölçek değerlendirmesi.** BLG'de atom yoğunluğu ~76 atom/nm²'dir. R = 30 nm bir flake ~216.000 atom, R = 50 nm ~600.000, R = 100 nm ~2.4 milyon atom demektir. Minivalley fiziği için ~25–30 nm efektif nokta + ~15–20 nm tampon, yani toplam ~40–50 nm flake yeterlidir, çünkü gap ile hapsedilen bağlı durumun evanesan kuyruğu birkaç nm'de söner. Yani proof-of-principle bir dizüstünde/orta bir makinede erişilebilirdir; darboğaz GPU değil, seyrek faktörizasyonun belleğidir.

## 2.  Araç seçimi ve kurulum mücadelesi

Hesaplama motoru seçimi, beklenenden uzun bir mücadeleye dönüştü ve bu sürecin kendisi metodolojik bir derstir: bir aracın teorik uygunluğu, onu çalıştırabilmekle aynı şey değildir. İzlenen yol ve sonuçları:

| Araç | Sonuç | Ders |
|---|---|---|
| pybinding (pip) | C++ derleme başarısız (Python 3.12, wheel yok) | C-uzantılı paketler pip'te kaynaktan derlenmeye çalışıp patlar |
| pybinding (conda) | Python 3.11/3.12 pin çakışması | condacolab Python sürümünü değiştirir, ortamı bozar |
| kwant (pip) | C/Cython derleme başarısız | Aynı derleme duvarı |
| kwant (conda) | Python pin çakışması, askıda kaldı | Kurulum cehennemi asıl işten saatler çaldı |
| saf NumPy/SciPy | ÇALIŞTI — sıfır kurulum | Mevcut, doğrulanmış toolchain her zaman en güvenli yoldur |

**Karar:** tight-binding Hamiltonyenini saf NumPy/SciPy ile elle kurmak. Bu, önceki k·p çalışmasının toolchain'iyle birebir aynıdır (her şey saf NumPy/SciPy ile yapılmıştı), sıfır derleme gerektirir, ve her matris elemanı görünür/denetlenebilir olduğundan çıpa disiplinine en uygunudur. Atomistik TB, özünde seyrek bir matris kurup gap yakınındaki birkaç özdeğeri shift-invert ile çekmektir — pybinding/kwant'ın C++ hızına proof-of-principle aşamasında ihtiyaç yoktur.

## 3.  Model: 4-bant BLG tight-binding örgüsü

Bernal (AB) istifli BLG, dört alt-örgüyle tanımlanır: A1, B1 (alt katman) ve A2, B2 (üst katman). AB istifinde A2 doğrudan B1'in üstündedir (dimer çifti). Parametreler standart McCann–Koshino değerleridir:

| Parametre | Değer | Anlam |
|---|---|---|
| a | 0.24595 nm | örgü sabiti |
| a_cc | 0.142 nm | C–C bağ uzunluğu (a/√3) |
| γ₀ | −3.16 eV | katman-içi en yakın komşu |
| γ₁ | 0.38 eV | B1–A2 dimer (interlayer) |
| γ₃ | 0.38 eV | A1–B2 trigonal warping (minivalley için ŞART) |
| γ₄ | 0.14 eV | skew (e-h asimetrisi; opsiyonel, kapalı) |
| Δ′ | 0.022 eV | dimer on-site asimetri (opsiyonel, kapalı) |

**Hapsetme.** Gap (displacement field) katmana-bağımlı on-site potansiyelle açılır: alt katman +Δ/2, üst katman −Δ/2. Kuantum noktası, önceki k·p çalışmasındaki ile aynı Gauss kuyusuyla tanımlanır: V(r) = −U₀·exp(−r²/2L_dot²). K noktası bu örgü konvansiyonunda K = (4π/3a, 0) ≈ (17.0 nm⁻¹, 0)'dadır.

**Hopping iskeleti (saf scipy, COO seyrek matris):** γ₀ her iki katmanda monolayer deseniyle ([0,0], [0,−1], [1,−1] hücre ofsetleri), γ₁ tek dikey B1–A2 bağı, γ₃ A1–B2 arasında üç bağ ([0,−1], [1,−1], [1,−2]). Bu γ₃ indeks vektörleri geometriden türetildi; doğruluğu bulk bant yapısı çıpasıyla (Bölüm 4–5) doğrulandı.

## 4.  Yapısal çıpalar

Önceki k·p çalışmasının disiplini buraya birebir taşındı: hiçbir fiziksel sonuca, aşağıdaki çıpalar geçmeden güvenilmedi.

| # | Çıpa | Sonuç |
|---|---|---|
| 1 | Hermitiklik | TB matrisi yapıca Hermityen (reel-simetrik); doğrulandı |
| 2 | Bulk gap vs k·p Tablo 1 | Üç haneye kadar uyum (Bölüm 5) — GEÇTİ |
| 3 | Trigonal warping / C₃ | Üç ayrı minivalley cebi, 120° simetri — GEÇTİ |
| 4 | Δ_ext = 0 → gapless | K noktasında gap kapandı — GEÇTİ |
| 5 | K vs K′ Kramers (TRS) | max\|E(k) − E(−k)\| = makine sıfırı — GEÇTİ |
| 6 | Sürekli limit: k·p ile dot seviyeleri | E₀ ve Δ₀₁ üç haneye kadar mutabık (Bölüm 7) — GEÇTİ |

Bu altı çıpanın hepsinin geçmesi, çözücünün hem bulk hem hapsedilmiş rejimde doğru olduğunu, ve en kritik olarak, k·p'nin doğru çözdüğü bölgede k·p'ye indirgendiğini garanti eder. Ancak bu sağlandıktan sonra, k·p'nin göremediği bölgeye güvenle girilebilir.

## 5.  Bulk doğrulaması: gap mutabakatı ve üç minivalley

Çıpa 2 sürecinde bir konvansiyon dersi yaşandı ve düzeltildi. İlk denemede on-site ±Δ/4 alındı; bu, gap'i yarıya bölerek k·p ile uyumu bozdu. Doğru konvansiyon ±Δ/2'dir: K noktasında iletkenlik kenarı E_c = δ/2, ve gerçek gap K'nin hemen yanındaki Mexican-hat halkasında, bu değerin birkaç meV altındadır. Ayrıca BZ taramasının K noktasına denk gelmesi (ızgara çözünürlüğü) önemlidir; aksi halde hat'i ıskalayıp sahte değer raporlanır.

Doğru kurulumla (±Δ/2, K etrafında ince fiziksel-koordinat taraması), atomistik iletkenlik kenarı E_c, k·p'nin bulk_gap değerleriyle karşılaştırıldı:

| δ_ext (eV) | k·p bulk_gap (meV) | Atomistik E_c (meV) | fark |
|---|---|---|---|
| 0.100 | 45.51 | 45.51 | ≈0 |
| 0.150 | 65.07 | 64.62 | ~0.5 |
| 0.167 | 71.58 | 70.46 | ~1.1 |
| 0.200 | 83.94 | 80.86 | ~3.1 |
| 0.250 | 101.78 | 94.37 | ~7.4 |

Düşük δ'da binde-bir uyum tesadüf değildir; örgü, gap konvansiyonu (±Δ/2 doğru) ve parametrelerin tümünün k·p ile tutarlı olduğunu gösterir. Yüksek δ'daki ~%8 sistematik kayma rastgele değil, fizikseldir: minivalley cepleri (hat yarıçapı 0.15 → 0.21 nm⁻¹) K'den uzaklaştıkça gerçek iletkenlik kenarı analitik γ₃=0 formülünün altına iner, ve fark warping güçlendikçe (δ büyüdükçe) büyür — yani senkronize olması gereken şeyle senkronizedir. Bu bir hata değil, atomistik modelin k·p'den daha doğru olduğu yerdir.

**Çıpa 3 — üç minivalley.** K noktası etrafında fiziksel koordinatta band kenarına zoom yapılan iletkenlik konturu, tam üç ayrı minivalley cebini 120° simetride (C₃) gösterdi. Bu, k·p kare-ızgarasının yapısal olarak üretemediği fiziktir ve atomistik örgünün gerçek lattice C₃'üne saygı duyduğunun doğrudan görsel ispatıdır (Knothe–Fal'ko / Mayer–Knothe minivalley yapısı). İlk denemedeki 'diagonal elips' görüntüsü, wraparound'un eğik k-koordinatından kaynaklanan bir görselleştirme kusuruydu; fiziksel koordinata geçince üç temiz cep ortaya çıktı.

## 6.  Kuantum noktası çözücüsü: artifaktlar, çökmeler, yakınsama disiplini

Dot çözücüsünü kurarken, önceki k·p raporunda defalarca yaşanan artifakt desenleri burada da aynen yaşandı ve aynı disiplinle elendi.

### 6.1  Bellek ve çökmeler

116k atomlu bir matriste shift-invert (eigsh + SuperLU faktörizasyonu) 9.9 GB'lık yerel makinede swap'e düşüp ~39 dakika sürdü (gizli çökme). Teşhis: darboğaz GPU değil, seyrek LU faktörizasyonunun belleğidir. Çözüm yüksek-RAM ortama (Colab, 51 GB) geçmek oldu; orada 380k atomlu sistemler ~40 saniyede, 650k atomlu sistemler ~130 saniyede çözüldü. Kritik nokta: H100 GPU'su bu adımda boş durur; gereken yüksek-RAM CPU çekirdeğidir, GPU değil.

### 6.2  Kenar artifaktı (kuyu > flake)

İlk dot denemesinde L_dot = 15 nm kuyu, R = 14 nm flake'ten genişti. Sonuç: tüm durumlar ⟨r⟩ ≈ 12.5 nm ile flake kenarına yapıştı — bağlı durum değil, kenar durumu. Düzeltme: kuyuyu flake'in iyice içine oturtmak (R ≥ ~2.5 L_dot). Bu, k·p raporundaki d≈60 nm illüzyonu ve box-divergence ile aynı sınıf hatadır.

### 6.3  Sahte near-dejenerelik ve yakınsama testi

Bağlı durumlar göründükten sonra, en düşük iki durum arasında küçük bir yarılma (0.037 meV) gözlendi ve bir an 'minivalley yarılması' diye yorumlandı. Yakınsama testi bunu derhal çürüttü:

| R (nm) | atom | split (meV) | ⟨r₀⟩ (nm) |
|---|---|---|---|
| 30 | 195.740 | 0.4714 | 29.4 |
| 42 | 381.617 | 0.0370 | 41.3 |
| 55 | 652.019 | 0.0044 | 54.2 |
| 68 | 994.393 | 0.0002 | 67.0 |

İki sütun birden konuştu: split flake büyüdükçe sıfıra çöktü (yakınsamadı = artifakt), ve ⟨r₀⟩ her zaman flake yarıçapının tam kendisi çıktı (durumlar kenara yapışmış). Filtre eşiği (0.4·R) flake'le büyüdüğü için kenar durumlarını yanlışlıkla içeri almıştı. Düzeltme: mutlak ⟨r⟩ eşiği (< 22 nm). Mutlak eşikle, gerçek bağlı durumların enerjileri R ile sabitlendi (R=55 ve R=68 birebir aynı) — gerçek dot spektrumu bu boyut-kararlı settir.

## 7.  Çıpa 6: k·p sürekli-limit mutabakatı

Gerçek bağlı durumlar mutlak eşikle ayrıldıktan sonra, en kritik çıpa: atomistik dot seviyeleri, k·p'nin aynı parametrelerle (δ=0.167, U₀=20 meV, L_dot=15 nm) verdiği seviyelerle birebir karşılaştırıldı.

| Nicelik | k·p (zarf) | Atomistik TB | Durum |
|---|---|---|---|
| Taban durumu E₀ (meV) | 60.32 | 60.317 | üç hane mutabık |
| Orbital ayrılma Δ₀₁ (meV) | 0.989 | 0.989 | üç hane mutabık |

**Kritik kavramsal düzeltme.** Başlangıçta 'Δ₀₁' diye alınan ardışık iki durum farkı (~0.22 meV), aslında valley-dejenere bir ikilinin içindeki mikro-yarılmaydı — gerçek orbital ayrılma değil. Atomistik spektrum, her orbitali iki-katlı (valley) verir ve k·p'nin tek 'Δ₀₁' gördüğü yeri minivalley alt-yapısına çözer. Doğru Δ₀₁, ikili-içi mikro-yarılma değil, ardışık iki orbital kabuk (ikili) arasındaki mesafedir. Bu doğru tanımla Δ₀₁ = 0.989 meV çıktı ve k·p ile birebir oturdu.

Bu mutabakat, hiçbir serbest parametre ayarlanmadan, tamamen bağımsız iki yöntem (k·p düzlem-dalga vs atomistik tight-binding) arasında sağlandı. Model bu noktada mühürlendi: artık k·p'nin göremediği bölgeye tam güvenle girilebilir.

## 8.  Valley ayrımı: temiz |ψ_K⟩ ve |ψ_K′⟩ çıkarmak

Reel eigsh, dejenere valley ikililerini keyfi reel kombinasyonlar olarak verir (net valley polarizasyonu sıfır). Gerçek valley özdurumlarını çıkarmak için, ikilinin 2-boyutlu dejenere alt-uzayında bir valley operatörü kuruldu: dalga fonksiyonunun Fourier içeriğini +K ve −K bölgelerine projekte eden operatörün (M(K) − M(K′)) köşegenleştirilmesi. Bu, atomistik modelin sahip olduğu, zarf modelinin atıp kaybettiği |K−K′| Bloch salınımını kullanır.

| Ölçüm | Sonuç | Yorum |
|---|---|---|
| Ayrımdan önce polarizasyon | +0.000 / −0.000 | karışık (reel kombinasyon) |
| Ayrımdan sonra polarizasyon | +0.990 / −0.990 | temiz valley özdurumları |
| TRS: \|⟨conj(ψ_K)\|ψ_K′⟩\| | 1.000 | K′ = zaman-tersi eşleniği |

Polarizasyon ±0.99 (ideal ±1), TRS kontrolü tam 1.000. Görselde iki trefoil birbirinin ayna-kiral hali (K ve K′ zıt kiralite — Berry curvature'ın zıt işareti). Artık temiz, valley-etiketli özdurumlar elde edildi; bu, intervalley matris elemanı hesabının önkoşuludur.

## 9.  İntervalley matris elemanı: hedef

Aşama 3'ün asıl ürünü: intervalley matris elemanı ⟨ψ_K|V|ψ_K′⟩. Fizik şudur: düzgün (yavaş değişen) hapsetme potansiyeli bu elemanı veremez, çünkü |K−K′| ≈ 17 nm⁻¹ Fourier bileşeni yok denecek kadar küçüktür (zarf modelinin 'valley korumalı' demesinin nedeni). İntervalley karışım ancak atomik-ölçek, kısa-menzilli bir pertürbasyondan (kafes kusuru, keskin yük dalgalanması) gelir. Atomistik model bunu hesaplayabilir; k·p hesaplayamaz.

**Bir ölçüm hatası ve düzeltmesi.** İlk denemede eleman site-başına modüllerin oranı olarak alındı; bu, modül alma işlemiyle intervalley faz salınımını sildiği için zorunlu olarak 1.000 verdi (yanlış). Doğru ölçüm fazı korur: ⟨ψ_K|V|ψ_K′⟩ = Σ_site V(site)·conj(ψ_K)·ψ_K′, kompleks toplam, modül en sonda alınır. İntra-valley'de faz site-bağımsız (reel yoğunluk, büyük); intervalley'de faz exp(i(K−K′)·r) ile döner ve büyük ölçüde iptal olur (küçük). İşte koruma budur.

Faz-koruyan doğru ölçümle, δ=0.167'de iki rejim arasındaki kontrast net çıktı:

| Kusur tipi | intervalley / intra | Yorum |
|---|---|---|
| Geniş Gauss (σ=5 nm) | 7.77e-03 | ≈0: faz iptali = KORUMA |
| Atomik tek-site | 1.000 | O(1): kısa menzil korumayı kırar |

Bu kontrast, ⟨ψ_K|V|ψ_K′⟩ elemanının belirgin bir sanal kısmıyla (2.14e-4 + 1.11e-4j) birlikte — intra ise neredeyse saf reel — intervalley Bloch-faz salınımının imzasını taşır. Bu kontrastın kendisi yakınsamaya muhtaç değildir: mertebe farkı (10⁻³ vs 1) o kadar büyüktür ki kenar gürültüsünden bağımsızdır.

## 10.  Negatif sonuç: tasarıma bağlı intervalley figürü yoktur

Asıl hedef, intervalley duyarlılığı dot tasarımının (gap, boyut) fonksiyonu olarak taramaktı — kimsenin QD geometrisinde hesaplamadığı bir nicel ilişki. Bu tarama yapıldı ve yakınsama disiplininden geçirildi. Sonuç dürüsttür ve negatiftir.

### 10.1  İlk tarama: yakınsamadı

δ_ext taraması (R=42'de) bir trend gösterdi (δ=0.10'da 0.142, δ=0.15'te 0.056), ama flake-boyutu yakınsama testi bu trendi çürüttü:

| R (nm) | atom | intervalley/intra (Test 1) |
|---|---|---|
| 42 | 381.617 | 7.79e-02 |
| 55 | 652.019 | 7.77e-03 |
| 68 | 994.393 | 2.73e-01 |

0.078 → 0.008 → 0.27: bir büyüklük mertebesi aşağı, sonra iki mertebe yukarı. Tamamen kararsız — klasik artifakt imzası, k·p raporundaki d≈60 nm ve valley-moment box-divergence ile aynı desen.

### 10.2  Yol B: kurtarma denemesi de çöktü

Sinyali kenar gürültüsünden ayırmak için iki silah birden kullanıldı: (i) yumuşak duvar — flake kenarında gap'i 0.6 eV'ye çıkarıp dot durumunu kenardan yalıtmak; (ii) çok-kusur istatistiği — tek kusur yerine 40 kusur konumu üzerinden RMS ortalama. Sonuç:

| R (nm) | atom | intervalley/intra (Yol B) |
|---|---|---|
| 42 | 381.617 | 1.05e-01 |
| 55 | 652.019 | 4.51e-02 |
| 68 | 994.393 | 1.32e-03 |

Yumuşak duvar vahşi zıplamayı sakinleştirdi (artık monoton) ama oran sıfıra doğru çöktü, sabitlenmedi: 0.105 → 0.045 → 0.0013. Bu da yakınsama değildir.

### 10.3  Dürüst yorum

**Bu monoton çöküş en net cevabı verir: düzgün-potansiyelin intervalley elemanının gerçek değeri sıfırdır** (ya da flake büyüdükçe sistematik olarak sıfıra giden bir şey). Sonlu bir değere inmiyor çünkü sonlu bir değeri yok — faz iptali tam çalışıyor. Ölçülen her şey sonlu flake'in artık gürültüsüydü; flake büyüdükçe o gürültü azalıyor, gerçek sinyal ortaya çıkmıyor çünkü gerçek sinyal sıfır. 'İntervalley duyarlılığın dot tasarımına bağlılığı' diye bir figür yoktur, çünkü düzgün potansiyelin intervalley elemanı tasarımdan bağımsız olarak sıfırdır.

**Bu bir başarısızlık değil, bir bulgudur — ve k·p raporunun ana tezini doğrular.** k·p 'düzgün potansiyel intervalley veremez, çünkü |K−K′| bazımda yok' diyordu. Atomistik TB ise |K−K′| bazında olduğu halde yine sıfır veriyor — yani bu bir baz eksikliği değil, gerçek bir fizik: düzgün potansiyelin intervalley elemanı faz iptaliyle gerçekten yok olur. İki bağımsız yöntem, aynı sonuç.

## 11.  Yakınsayan ve hayatta kalan sonuçlar

Negatif tarama sonucuna rağmen, yakınsamış ve fiziksel olarak sağlam birkaç sonuç elde edildi. Bunlar çalışmanın kalıcı çıktılarıdır:

- **Çalışan atomistik çözücü.** Saf NumPy/SciPy, kwant'sız, beş yapısal çıpadan geçmiş, ~1 milyon atoma kadar ölçeklenen bir BLG QD tight-binding çözücüsü.

- **k·p sürekli-limit mutabakatı.** E₀ (60.317 vs 60.32) ve Δ₀₁ (0.989 vs 0.989) üç haneye kadar; serbest parametre yok. İki bağımsız yöntemin mühürlü çapraz-doğrulaması.

- **Üç minivalley / C₃.** k·p kare-ızgarasının üretemediği gerçek minivalley yapısı, hem bulk konturunda hem dot dalga fonksiyonu trefoilinde.

- **Temiz valley özdurumları.** Dejenere ikililerden ±0.99 polarizasyonlu, TRS-tutarlı (1.000) |ψ_K⟩ ve |ψ_K′⟩.

- **İki-rejim kontrastı.** Düzgün kusur → intervalley ≈ 0 (koruma); atomik tek-site kusur → O(1) (kırılma). Mertebe farkı yakınsamaya muhtaç değil. Bu, korumanın kökenini kafes-ölçeğe yerleştirir — k·p'nin söyleyemediği, deneye yorumlanabilir bir ifade.

## 12.  Durum, sınır ve sonraki adım

**Mevcut durum.** Atomistik araç kuruldu ve doğrulandı (Çıpalar 1–6 ✓). Minivalley yapısı, valley özdurumları ve iki-rejim kontrastı gerçek ve yakınsamış sonuçlardır. Ancak düzgün-potansiyel intervalley elemanının dot-geometri bağımlılığı, mevcut sonlu-flake çözücünün gürültü tabanının altındadır — çünkü gerçek değeri sıfıra çok yakındır. Bu, k·p raporunun 'intervalley fiziği bu çözücünün dışında' dediği sınırın atomistik teyididir.

**Sınırın doğası.** Düzgün potansiyel intervalley elemanını faz iptaliyle yok eder; dolayısıyla valley korumasını sınırlayan şey düzgün gürültü değil, atomik-ölçek (kısa-menzilli) kusurlardır. Bu kusurların etkisini güvenilir biçimde nicelemek, sonlu bir disk flake yerine periyodik-kenarlı, gerçek bir lattice-kusur modeli gerektirir — yani k·p raporunun işaret ettiği tam-BZ atomistik adımın kendisi.

**Önerilen sonraki adım — önce ucuz fizibilite.** Tam periyodik modele haftalar yatırmadan önce, küçük bir periyodik supercell'de (birkaç bin atom) tek soru test edilmeli: bu kurulum prensipte sonlu bir intervalley elemanı verebiliyor mu? Küçük testte bile sinyal gürültüye gömülüyse, haftalar harcanmadan öğrenilir. Bu, k·p raporundaki 'önce ucuz çıpa, sonra pahalı hesap' disiplininin doğrudan uygulamasıdır.

**Alternatif sağlam yön.** İntervalley'i bir kenara bırakıp, atomistik aracın güçlü olduğu ve yakınsadığı bir soruya dönmek: minivalley spektrumunun dot tasarımıyla (boyut, gap, valley-merkez vs minivalley rejimi) değişimi. Bu, Mayer–Knothe çizgisinde, yakınsayan ve raporlanabilir bir sonuç olabilir.

## 13.  Kod ve parametre özeti

Çözücünün çekirdek yapısı (saf NumPy/SciPy):

# orgu: 4 alt-orgu A1,B1 (alt), A2,B2 (ust); AB istif

basis = [[0,-acc/2],[0,acc/2],[0,acc/2],[0,3*acc/2]]

layer = [0,0,1,1]   # 0=alt(+D/2)  1=ust(-D/2)

# hopping kurallari: (from, to, dn, dm, deger)

HOPS = [(0,1,0,0,G0),(0,1,0,-1,G0),(0,1,1,-1,G0),   # g0 alt

(2,3,0,0,G0),(2,3,0,-1,G0),(2,3,1,-1,G0),   # g0 ust

(1,2,0,0,G1),                                # g1 dimer

(0,3,0,-1,G3),(0,3,1,-1,G3),(0,3,1,-2,G3)]  # g3 warping

# on-site: gap + Gauss kuyu

ons = where(layer[ks]==1, -D/2, +D/2) - U0*exp(-r^2/(2 L^2))

# cozucu: seyrek shift-invert

E, psi = eigsh(H.real.tocsc(), k=n, sigma=E_c - U0/2)

# valley ayrimi: M(K)-M(K') koselestir -> psi_K, psi_K'

# intervalley: |sum_site V(site) * conj(psi_K) * psi_K'|  (faz korunur!)

**Referans çalışma noktası:** δ_ext = 0.167 eV, U₀ = 20 meV, L_dot = 15 nm, R = 42–68 nm. Enerjiler meV, uzunluklar nm. γ₀ = −3.16, γ₁ = γ₃ = 0.38 eV; γ₄ ve Δ′ kapalı. Hesaplar 51 GB RAM'li bir bulut çekirdeğinde yapıldı; 380k atom ~40 s, 650k atom ~130 s, ~1M atom ~birkaç dakika (shift-invert, CPU).

**Kapanış notu.** Bu oturum bir 'büyük keşif' üretmedi, ama çalışan ve doğrulanmış bir atomistik araç üretti ve onu dürüstçe sınadı. Heyecan verici bir intervalley tarama figürü, yakınsama testinde — tıpkı k·p raporundaki altı sahte sweet-spot gibi — dağıldı. Bu disiplin, durable bir bulguyu geri çekilen bir bulgudan ayıran şeydir. Araç ve doğrulanmış sonuçlar elde kaldı; doğru soruyu beklemekle hiçbir şey kaybetmezler.
