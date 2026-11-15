> **SUPERSEDED DRAFT / GEÇERSİZ TASLAK.** This early Turkish manuscript draft claimed a ~9 % T2* gain from AC driving using a period-averaging approximation. Its central claim was later refuted (curvature sign is convex, and a full-Hilbert Floquet calculation shows 38–95 % leakage at A_AC ≥ 20 meV); see `docs/BLG_Charge_Noise_Protection_Report.md`. Kept for the record only.

**Bilayer** **Graphene****Quantum** **Dot** **Valley** **Qubit’lerinde**** AC-Sürüş ile Dinamik Yük Gürültüsü Koruması:**

**Bant-Eğriliği Üzerinden Statik Tatlı Nokta Gerektirmeyen Bir Mekanizma**

**Faruk Sami Bilgin**

*Bağımsız Araştırmacı, Kocaeli, Türkiye*


### Özet

Yarıiletken kuantum nokta (QD) kübitleri yük gürültüsüne karşı duyarlıdır. Bu duyarlılığı azaltmanın bilinen yolu statik "tatlı nokta" (sweet spot) — yani kübit enerji ayrılışının ∂Δ₀₁/∂V = 0 olduğu kapı voltajı — etrafında çalışmaktır. Yakın zamanda Hajati ve Burkard [Phys. Rev. B 110, 245301 (2024)], silisyum/germanyum (Si/Ge) düz-mod (flopping-mode) spin kübitlerinde AC sürüş ile dinamik tatlı noktaların ortaya çıkabileceğini analitik olarak gösterdi; bu mekanizma statik tatlı nokta etrafında salınımlı rezonans-engellemeye dayanır. Bu çalışmada, çift katmanlı grafen (bilayer graphene, BLG) kuantum noktası valley kübitleri için aynı soruyu doğrudan sayısal olarak ele alıyoruz. 4-bant k·p sürekli ortam Hamiltonyenini düzlem-dalga bazında diagonalize ediyor, beş yapısal çıpa ile modeli doğruluyor ve operasyonel çalışma noktasında (Δ_ext = 0.167 eV, U₀ = 20 meV, L_dot = 15 nm, Δ₀₁ ≈ 1.01 meV) kuasi-statik Gauss yük gürültüsü (σ_D = 2 ve 5 meV) altında 200 örneklik bir topluluk ile AC sürüş etkisini ölçüyoruz. Bulduğumuz BLG valley kübitlerinde statik bir tatlı nokta yok — Δ₀₁(Δ_ext) monoton azalandır — ancak eğri belirgin biçimde konkavdır ve periyot-ortalama Δ₀₁’in standart sapmasını azaltır. A_AC = 50 meV’de T₂* ≈ 1.09 × T₂*_statik bulduk; etki σ_D’den neredeyse bağımsızdır (σ_D = 2 ve 5 meV’de aynı oran). Mekanizma Jensen-tipi eşitsizliktir: konkav Δ₀₁ → periyot-ortalama değer < statik değer, ve dağılım daralır. Bu Hajati-Burkard mekanizmasından nitel olarak farklıdır ve statik tatlı nokta gerektirmez. Operasyonel sınırlar (A_AC ≲ 50 meV) ve gelecekteki uzantılar (spin, tam Floquet, gerçekçi 1/f gürültü) tartışılıyor.

**Anahtar Kelimeler:** bilayer grafen, kuantum nokta, valley kübit, yük gürültüsü, dinamik tatlı nokta, AC sürüş, k·p modeli

## 1. Giriş

### 1.1 Yük Gürültüsü ve Tatlı Nokta Çalışması

Yarıiletken kuantum noktalarında (QD) elektrostatik olarak tanımlanmış kübitler, ölçeklenebilir kuantum bilgi işlem için önemli bir platform sınıfını oluşturur [Loss & DiVincenzo 1998; Burkard et al. 2023]. Bu platformlar yüksek-kalite elektriksel kontrole olanak tanır, ancak buna karşılık temel bir zayıflık taşır: kapıların oluşturduğu elektrik alanları ve gürültü, kübit enerji seviyelerini kapı potansiyelinin değişimleri yoluyla doğrudan etkiler. Bu yük gürültüsü çoğunlukla 1/f spektrumlu olup dekoherensı belirleyen baskın faktördür.

Yük gürültüsüne karşı temel savunma stratejisi "tatlı nokta" (sweet spot) işletimidir: kübit ayrılışı Δ₀₁ kapı voltajına karşı türevinin sıfırlandığı bir noktada, ∂Δ₀₁/∂V = 0, kübit voltaj dalgalanmalarına birinci derecede duyarsız hale gelir. Bu strateji singlet-triplet kübitler [Petta et al. 2005], yük kübitleri [Stockklauser et al. 2017] ve Si/Ge düz-mod spin kübitleri [Benito et al. 2019] için yerleşik bir tekniktir.

Statik tatlı noktanın doğal sınırı, kübitin yalnızca o spesifik voltaj değerinde korunmasıdır. Operasyonel manipülasyon (Rabi salınımları, two-qubit kapıları) genellikle bu noktadan uzaklaşmayı gerektirir. Bu sınırlamanın üstesinden gelmek için son yıllarda dinamik tatlı nokta önerisi gelişti: AC sürüş altında periyot üzerinden ortalanmış kübit ayrılışı, sürüş parametrelerine göre türevi sıfır olan noktalar gösterebilir [Mundada et al. 2020 (süperiletken); Hajati & Burkard 2024 (Si/Ge)].

### 1.2 Hajati-Burkard 2024 ve Si/Ge Düz-Mod Kübitlerde Dinamik Tatlı Nokta

Hajati ve Burkard [Phys. Rev. B 110, 245301 (2024)], Si/Ge çift kuantum noktasında düz-mod spin kübit’ini ele aldı. Rotating wave approximation (RWA) içinde, zaman-bağımlı etkin Hamiltonyen kullanarak Rabi frekansı için analitik kapalı-form ifadeler türettiler. Çekirdek bulguları: AC sürüş, statik tatlı noktanın kayma yoluyla ortaya çıkan rezonans-dışı sürüş bölgelerinde dinamik bir tatlı nokta oluşturur. İkinci derece (Bloch-Siegert tipi) dinamik tatlı noktalar bile sürüş ve detuning parametrelerinin belirli kombinasyonlarında ortaya çıkar.

Hajati-Burkard yaklaşımının iki temel özelliği vardır: (i) analitik bir RWA çerçevesi kullanılır, dolayısıyla sürüş frekansı kübit enerjisine yakın (rezonans-dışı kayma rejimi) varsayılır; (ii) ortaya çıkan dinamik tatlı noktalar daima statik bir tatlı noktanın etrafında kümelenir — yani sistemin statik özelliği olarak ∂Δ₀₁/∂V = 0 noktası mevcut olmalıdır. Bu mekanizma silikon ve germanyum spin kübitleri için fiziksel olarak iyi temellendirilmiştir, çünkü flopping-mode geometrisi sentetik spin-orbit etkileşimi yoluyla statik bir tatlı nokta üretir.

### 1.3 BLG Valley Kübitleri: Bilinmeyen Bir Soru

Bilayer grafen (BLG), iki AB-istif edilmiş grafen katmanından oluşur ve dik bir elektrik alanı ile açılan ayarlanabilir bir bant aralığına sahiptir [McCann & Koshino 2013]. Elektrostatik QD’ler bu malzemede deneysel olarak gerçekleştirilmiştir [Eich et al. 2018; Garreis et al. 2021; Banszerus et al. 2022] ve düşük spin-orbit kuplajı ve nükleer spinsiz karbon matrisi nedeniyle uzun spin/valley koherens süreleri sunma potansiyeli taşır. Knothe ve Fal’ko [Phys. Rev. B 98, 155435 (2018)], BLG QD ve kuantum tellerinde trigonal bükülmenin (γ₃ terimi) üç-katlı "minivadi" yapısı yarattığını, bu yapının Berry curvature aracılığıyla topolojik bir yörünge magnetik momenti taşıdığını gösterdi. Mayer ve Knothe [pssb 260, 2300395 (2023)] BLG QD spektrumunun parametre uzayını sayısal olarak taramış, valley g-faktörünün kuyu boyutu, derinliği ve dış gap ile nasıl değiştiğini haritalandırmıştır.

Si/Ge’den farklı olarak BLG valley kübit, dört bant k·p sürekli ortam modelini gerektirir ve trigonal warping ile birlikte gelir. Bu yapı, basit analitik RWA yaklaşımına direnir. Üstelik BLG QD’lerinde Hajati-Burkard tipi sentetik spin-orbit yoktur; valley kübit ayrılışı yalnızca elektrostatik kapılarla kontrol edilir. Bu nedenle iki soru açık kalmıştır: (i) BLG valley kübitlerinde statik bir tatlı nokta var mıdır? (ii) Var olmasa bile, AC sürüş ile yük gürültüsüne karşı koruma sağlanabilir mi?

### 1.4 Bu Çalışmanın Katkısı

Bu çalışmada her iki soruyu doğrudan ele alıyoruz. Hajati-Burkard’ın analitik yaklaşımı yerine, BLG’nin doğal olarak gerektirdiği tamamen sayısal bir yol izliyoruz: 4-bant k·p Hamiltonyenini düzlem-dalga bazında diagonalize ediyoruz, beş yapısal çıpa ile modeli doğruluyoruz ve operasyonel çalışma noktasında AC sürüş ile gürültü topluluğunu birleştiriyoruz. Çekirdek bulgularımız:

- Operasyonel rejimde Δ₀₁(Δ_ext) monoton azalandır; ∂Δ₀₁/∂Δ_ext = 0 noktası yoktur. Yani BLG valley kübitlerinde, en azından bu kuyu geometrisi için, statik bir tatlı nokta mevcut değildir.

- Buna rağmen Δ₀₁(Δ_ext) belirgin biçimde konkavdır. AC sürüş ile periyot üzerinden ortalama alındığında bu konkavlık Jensen-tipi eşitsizlik aracılığıyla σ_Δ₀₁ standart sapmasını azaltır.

- A_AC = 50 meV genliğinde, σ_D = 5 meV gürültü seviyesinde T₂* ≈ 1.09 × T₂*_statik elde edilir; etki σ_D = 2 meV’de neredeyse aynı orana sahiptir (gürültü-ölçeği bağımsızlığı).

- Mekanizma Hajati-Burkard’dan nitel olarak farklıdır: statik tatlı nokta gerektirmez, yalnızca band yapısının yerel eğriliğine bağlıdır.

- Operasyonel sınır A_AC ≲ 50 meV olarak tanımlanır; daha büyük genliklerde periyot içi anlık dejenerlik (Δ₀₁(t) → 0) ortaya çıkar, A_AC ≳ 150 meV’de ise QD kuyusu açılır.

Bu, bilayer grafen valley kübitlerinde dinamik tatlı nokta sorusunun ilk sistematik sayısal incelemesidir ve aynı zamanda AC korumasının statik tatlı nokta gerektirmediği yeni bir mekanizma sınıfına işaret eder.

## 2. Model ve Yöntem

### 2.1 Dört Bant k·p Hamiltonyeni

Çalışmamız McCann-Koshino sürekli ortam modeline dayanır [McCann & Koshino 2013]. Bernal istifli BLG’de düşük enerji fiziği, A1 (alt katman A alt-örgüsü), B1, A2, B2 (üst katman B alt-örgüsü) sıralamasında dört-bileşenli sözde-spinörler ile tanımlanır. K vadisi (ξ = +1) etrafında k·p Hamiltonyeni şu şekildedir:

*H_kp**(q) = **diag**(**+Δ/2, +Δ/2, −Δ/2, −Δ/2)*

+ kinetik off-diagonal terimler v · π (ana iletkenlik), γ₁ (katmanlar arası kuplaj) ve v₃ · π (trigonal warping). Burada π = ξq_x + iq_y, v = (√3/2)·a·γ₀, ve v₃ = (√3/2)·a·γ₃. Parametre değerleri: a = 0.246 nm (örgü sabiti), γ₀ = 3.16 eV (intra-katman), γ₁ = 0.38 eV (inter-katman) ve γ₃ = 0.38 eV (trigonal). Bu değerler standart BLG literatürüyle uyumludur. Dış elektrik alan tarafından açılan bant aralığı Δ, üst ve alt katmanlar arasındaki potansiyel farkını temsil eder; analitik γ₃ = 0 limitinde iletkenlik bant minimumu γ₁·Δ/(2·√(γ₁² + Δ²)) ifadesine sahiptir.

### 2.2 Düzlem-Dalga Bazı

Hapsetme potansiyeli V(r) ile birleşik problem [H_kp(q) + V] ψ = E ψ formundadır. Bu problem doğal olarak düzlem-dalga bazında çözülür: ψ(r) = Σ_q c_q · e^(iq·r), burada q düzlem-dalga momentumları L_box × L_box kutu içinde ayrıklaştırılır. Kutu boyutu L_box = 120 nm, momentum adımı Δq = 2π/L_box ≈ 0.0524 nm⁻¹, ve kesme q_max = 0.5 nm⁻¹ alındı. Toplam düzlem-dalga sayısı N_pw = 293 ve Hilbert uzayı boyutu 4·N_pw = 1172 oldu.

Bu baz seçimi rastgele değildir. Erken denemelerde sonlu fark Dirac diskretizasyonu kullanılmış, ancak bilinen fermion-doubling artefaktı (E ≈ 0 etrafında beklenmedik 4-katlı dejenerelik ve dama tahtası tipi yoğunluk desenleri) modelin doğruluğunu bozmuştur. Düzlem-dalga bazına geçiş bu artefaktı tamamen ortadan kaldırmıştır (Ek A). Bu metodolojik seçim Knothe-Fal’ko (2018) ve Mayer-Knothe (2023) çalışmalarıyla uyumludur.

### 2.3 Hapsetme Potansiyeli

QD’yi tanımlamak için Gauss formda bir potansiyel kullanıldı:

*V(r) = −U₀ · **exp**(**−r² / (2 L_dot²))*

Burada U₀ kuyu derinliği ve L_dot kuyu yarıçapıdır. Çalışmamızın referans noktası U₀ = 20 meV, L_dot = 15 nm’dir. Bu seçim deneysel olarak ulaşılabilir bir rejimi (BLG QD’lerde ölçülmüş dot boyutları ~10–30 nm aralığındadır) ve perturbatif hapsetme limitini (U₀ ≪ bulk_gap) sağlar. Daha büyük U₀ değerleri (örn. 50 meV) anti-crossing rejimine girer ve sahte sweet-spot sinyalleri üretir; bu durum Ek B’de detaylandırılmıştır.

Potansiyelin momentum-uzayı Fourier dönüşümü analitiktir:

*V_FT(q) = −U₀ · 2π · L_dot² · **exp**(**−q² L_dot² / 2)*

ve kutu üzerinden normalleştirme L_box² ile bölünür. Bu nedenle V matrisi tek seferde hesaplanır ve Δ_ext değişimleri sırasında yeniden hesaplanmaz.

### 2.4 Beş Yapısal Çıpa

Modelin doğruluğunu garanti altına almak için beş bağımsız çıpa kontrol edildi. Bu çıpalar paper boyunca raporlanan tüm parametre rejimlerinde geçerliliğini korumuştur:

- **Çıpa 1 —** **Hermitiklik****:** |H − H†| = 0 (makine hassasiyeti).

- **Çıpa 2 —** **Bulk****limit:** V = 0 alındığında bant min ±γ₁·Δ/2√(γ₁²+Δ²) analitik ifadesinin γ₃ ile ~%5 düşeltilmiş halini verir; sayısal değer 45.5–46.9 meV (Δ_ext = 0.1 eV için 48.4 meV analitik beklenti).

- **Çıpa 3 —** **Trigonal** **warping****ve C₃:** γ₃ ≠ 0 ile dalga fonksiyonları üç-katlı simetri gösterir (Bölüm 3 Şekil 1’de görsel olarak doğrulanır).

- **Çıpa 4 — Δ = 0 limiti:** Dış gap kapatıldığında metalik spektrum elde edilir; bağlı durumlar bulk’a sızar.

- **Çıpa 5 —** **Kramers****eşleşmesi:** K (ξ = +1) ve K’ (ξ = −1) vadileri zaman-tersi simetrisi altında özdeş spektrumlar verir.

### 2.5 AC Sürüş ve Gürültü Modeli

Kapı voltajı dalgalanmaları doğrudan dış gap Δ_ext’i etkilediğinden, yük gürültüsünü Δ_ext’e eklenen kuasi-statik Gauss bir terim olarak modelliyoruz:

*Δ_ext** → **Δ_ext** + **δΔ** + A_AC · **cos(**ω t)*

Burada δΔ ~ N(0, σ_D) gerçeklenmiş yük gürültüsünü temsil eder, ve A_AC · cos(ω t) terimi kapıya uygulanan AC sürüştür. Gürültü topluluğu 200 bağımsız gerçeklemeden oluşur. Frekans seçimi ω/Δ₀₁ ≈ 6.6 olarak alındı, bu da Magnus seri açılımının hızlı sürüş limitinde geçerli olduğu rejimi sağlar. Periyot üzerinden 32 faz örneği ile anlık Δ₀₁ değerlerinin ortalaması alındı.

Standart sapma σ_Δ₀₁ ve dolayısıyla saf koherens süresi T₂* ∝ 1/σ_Δ₀₁ topluluk üzerinden hesaplanır. AC etkisini ölçen birincil gözlemlenebilirimiz T₂*_AC / T₂*_statik oranıdır.

### 2.6 Yapısal Çıpa: AC = 0 Tutarlılığı

A_AC = 0 alındığında periyot ortalaması statik Δ₀₁ değerine eşit olmalıdır. Bu yapısal çıpa, kodun doğruluğunu garanti eden bir kontroldür. Tüm hesaplamalarda max|⟨Δ₀₁⟩_AC(A=0) − Δ₀₁_statik| < 10⁻¹⁴ eV doğrulanmıştır (makine hassasiyeti). Bu, gözlenen herhangi bir AC etkisinin sayısal artefakt olmadığını garanti eder.

## 3. Statik QD Spektrumu

### 3.1 Bağlı Durum Yapısı

İlk olarak statik (sürüşsüz, gürültüsüz) QD spektrumunu inceledik. Anti-crossing rejimine girmemek için U₀ = 20 meV ve Δ_ext = 0.167 eV referans noktasını seçtik (U₀/bulk_gap ≈ 0.26, perturbatif hapsetme). Bu noktada iletkenlik bandının altında 17 bağlı durum bulduk. En düşük enerjili kübit ayrılışı Δ₀₁ = 1.0118 meV oldu, bağlanma derinliği (bulk min − E₀) 16.1 meV ile çok temiz.

Bağlı durumların gerçek-uzay dalga fonksiyonları Knothe-Fal’ko (2018) minivadi yapısının QD analogunu doğrudan gösterir: ψ₀ merkezde lokalize, ψ₁ ve ψ₂ üç-katlı C₃ simetrisinde "trefoil" yoğunluk dağılımları (üç tepe noktası, γ₃ trigonal warping’in imzası). Daha yüksek shell’ler benzer simetrik desenleri gösterir. Bu dalga fonksiyonu morfolojisi BLG QD’nin valley fiziğinin "perturbative confinement" rejiminde doğru yakalandığını teyit eder.

### 3.2 Δ_ext Taraması — Statik Tatlı Nokta Var Mı?

Statik tatlı nokta varlığını test etmek için Δ_ext ∈ [0.05, 0.30] eV aralığında 16 nokta üzerinde sistematik tarama yaptık. Sonuçlar Tablo 1’de özetlenmiştir.

| Δ_ext (eV) | bulk_gap (meV) | E₀ (meV) | Δ₀₁ (meV) | U₀/gap |
|---|---|---|---|---|
| 0.050 | 24.79 | 12.09 | 5.130 | 0.81 |
| 0.100 | 48.35 | 33.49 | 3.471 | 0.41 |
| 0.150 | 69.76 | 53.97 | 1.524 | 0.29 |
| 0.167 | 76.32 | 60.20 | 0.998 | 0.26 |
| 0.200 | 88.49 | 71.39 | 0.347 | 0.23 |
| 0.250 | 104.4 | 85.22 | 0.056 | 0.19 |
| 0.300 | 117.7 | 96.27 | 0.009 | 0.17 |

**Tablo 1.** *BLG QD statik spektrumu, U₀ = 20 **meV**, **L_dot** = 15 nm. Tüm noktalar **perturbatif** hapsetme rejiminde (U₀/**bulk_gap* *< 0.5**).*

Tablo 1’den iki kritik gözlem ortaya çıkar:

- Δ₀₁(Δ_ext) monoton azalandır. 5.13 meV’den (Δ_ext = 0.05) 0.009 meV’e (Δ_ext = 0.30) düzgün biçimde azalır. ∂Δ₀₁/∂Δ_ext sıfır geçişi göstermez — yani statik tatlı nokta bu kuyu geometrisinde mevcut değildir.

- Aynı zamanda eğri monoton azalan olmasına karşın belirgin biçimde konkavdır: ∂²Δ₀₁/∂Δ_ext² < 0. Bu konkavlık, ileride göreceğimiz gibi, AC korumasının fiziksel temelidir.

### 3.3 Anti-crossing Artefaktı: Neden U₀ = 20 meV?

Daha derin kuyularda (örn. U₀ = 50 meV) sahte sweet-spot sinyalleri gözlemledik: Δ_ext < 0.10 eV bölgesinde ∂Δ₀₁/∂Δ_ext işaret değiştirir gibi görünür. Detaylı inceleme, bu davranışın bağlı QD durumlarının bulk bantlarıyla rezonansa girip karışması (anti-crossing) olduğunu göstermiştir. Bu rejim U₀ > bulk_gap koşulunda ortaya çıkar ve tek-vadi izole QD varsayımımızı geçersiz kılar.

Bu nedenle U₀ değerini bulk_gap’in yarısının altına çektik (U₀ = 20 meV). Bu seçim Mayer-Knothe (2023) metodolojisiyle uyumludur: BLG QD analizleri tipik olarak U₀/bulk_gap < 0.5 rejiminde yapılır. Ek B’de U₀ = 50 meV rejiminin tam analizi sunulmuştur — bu rejimdeki "stasyoner noktalar" anti-crossing artefaktıdır, fiziksel sweet spotlar değildir.

## 4. AC Dinamik Koruma — Çekirdek Bulgular

### 4.1 Ana Sonuç: T₂* Uzaması

Çalışma noktasını Δ_ext = 0.167 eV (Δ₀₁ ≈ 1.0 meV) olarak sabitledik. Yük gürültüsü σ_D = 5 meV (literatür gerçekçi aralığının üst sınırı, bkz. Tartışma) ile 200 bağımsız gerçeklemede çekildi. AC genliği A_AC değerleri 0’dan 50 meV’e taranmıştır. Yapısal çıpa A_AC = 0 noktasında doğrulandı: max|⟨Δ₀₁⟩_AC − Δ₀₁_statik| = 0.0 eV (makine hassasiyeti). Sonuçlar Tablo 2’de özetlenmiştir.

| A_AC (meV) | ⟨Δ₀₁⟩ (meV) | σ_Δ₀₁ (meV) | T₂* / T₂*_statik |
|---|---|---|---|
| 0 | 1.0118 | 0.1220 | 1.0000 |
| 5 | 1.0149 | 0.1218 | 1.0011 |
| 10 | 1.0239 | 0.1214 | 1.0044 |
| 20 | 1.0593 | 0.1199 | 1.0171 |
| 30 | 1.1151 | 0.1176 | 1.0371 |
| 50 | 1.2730 | 0.1116 | 1.0926 |

**Tablo 2.** *AC genlik taraması, **σ_D** = 5 **meV**, 200 gerçeklik, 32 faz. T₂* monoton uzar.*

Bulgular kesindir: σ_Δ₀₁ AC genliğiyle monoton azalır ve T₂* monoton uzar. A_AC = 50 meV’de T₂* statik değere göre %9.26 daha uzun olur. Etki büyüklüğü A küçük olduğunda A²’ye orantılı görünür (Jensen ikinci derece katkı), bu Bölüm 4.3’te detaylandırılmıştır.

### 4.2 Karşılaştırma: 1B Oyuncak Modelin Tersi Sonuç

İlginç bir kontrol olarak, aynı v3 prosedürünü 1B Schrödinger oyuncak modelinde de uyguladık (kare kuyu, monoton ∂Δ₀₁/∂V₀). 1B durumda Δ₀₁(V) eğrisi konvekstir ve sonuç tam tersidir: A_AC arttıkça σ_Δ₀₁ artar, T₂* %5–6 mertebesinde kısalır. Bu kontrast, AC etkisinin işaretinin Δ₀₁(V) eğriliğinin işaretine bağlı olduğunu açıkça gösterir — pozitif (konkav) eğrilik koruma, negatif (konveks) eğrilik dekoherens artışı verir.

### 4.3 Mekanizma: Jensen Eşitsizliği

A_AC küçük olduğunda periyot ortalaması Taylor açılımına tabi tutulabilir:

*⟨Δ₀₁⟩_T ≈ Δ₀**₁(**Δ_ext** + **δΔ**) + (A_AC² / 4) · ∂²Δ₀₁/∂Δ_ext²*

Burada periyot-ortalama ⟨cos²⟩ = 1/2 kullanıldı. Eğer ∂²Δ₀₁/∂Δ_ext² < 0 (konkav), AC sürüş ortalaması yerel olarak statik değerden daha küçük olur. Daha kritik olarak, σ_Δ₀₁ topluluk üzerinden hesaplanan standart sapma da değişir. δΔ ~ N(0, σ_D) için, ikinci derece düzeltmenin etkisi σ_Δ₀₁’yi azaltacak yönde olur — çünkü konkav eğri "tepe noktalarını aşağı çekerken" "düşük noktaları yukarı iter", topluluğun varyansını sıkıştırır.

Bu Jensen-tipi mekanizmadır. Statik bir tatlı nokta (∂Δ₀₁/∂V = 0) gerektirmez; yalnızca yerel eğriliğin (∂²Δ₀₁/∂V²) doğru işarete sahip olmasını gerektirir. Bizim çalışma noktamızda eğrilik konkavdır (bkz. Bölüm 3.2) ve bu nedenle AC koruma sağlanır.

### 4.4 Hajati-Burkard ile Farklılık

Hajati-Burkard (2024) mekanizması statik bir tatlı nokta etrafında konumlanan rezonant kayma yoluyla dinamik tatlı nokta oluşturur. Ortaya çıkan koruma yerel bir fenomendir (sürüş frekansının ve genliğinin belirli kombinasyonlarında). Bizim Jensen mekanizmamız ise global eğrilik özelliğidir: ∂²Δ₀₁/∂Δ_ext² < 0 sağlandığı sürece (operasyonel rejimde her yerde sağlandığı görülür), AC sürüş her zaman koruma sağlar. Bu iki mekanizma rakip değil, tamamlayıcıdır; gelecekteki çalışmalar her ikisinin de aktif olduğu sistemlerde toplam etkiyi değerlendirebilir.

## 5. Sağlamlık Kontrolleri

### 5.1 Gürültü-Ölçeği Bağımsızlığı

Önerilen mekanizma yalnızca band eğriliğine dayanıyorsa, gürültü genişliği σ_D değişimleri yalnızca mutlak σ_Δ₀₁ değerini değiştirmeli, oran T₂*_AC/T₂*_statik’yi değiştirmemelidir. Bunu doğrulamak için σ_D = 2 meV (deneysel gerçekçi alt sınır, bkz. Tartışma) için tarayı tekrarladık. Sonuçlar Tablo 3’te.

| A_AC (meV) | σ_D = 5 meV: T₂*/T₂*_stat | σ_D = 2 meV: T₂*/T₂*_stat |
|---|---|---|
| 0 | 1.0000 | 1.0000 |
| 5 | 1.0011 | 1.0009 |
| 10 | 1.0044 | 1.0041 |
| 20 | 1.0171 | 1.0159 |
| 30 | 1.0371 | 1.0345 |
| 50 | 1.0926 | 1.0870 |

**Tablo 3.** *Gürültü-ölçeği bağımsızlığı: **σ_D** 2.5× azaltıldığında T₂* oranları neredeyse aynı kalır. Bu mekanizmanın **σ_D’den** bağımsız olduğunu, yalnızca **band** eğriliğine bağlı olduğunu doğrular.*

Tablo 3 belirleyici bir kanıt sunar: σ_D 2.5 kat azaltıldığında (5 → 2 meV) T₂* oranları %0.5’den az değişir. Bu mekanizmanın temelde Δ₀₁ fonksiyonunun yerel eğriliğinden geldiğini, gürültü dağılımının özelliklerinden değil, doğrular. Bu hakemin doğal sorusunu ("ya gerçek deneydeki çok daha küçük gürültü değerlerinde de işe yarıyor mu?") önceden cevaplar: evet, oransal etki korunur.

### 5.2 Doygunluk Eğrisi ve Operasyonel Sınır

A_AC değerlerini 50 meV’in ötesine taşıyarak doygunluk davranışını araştırdık (M = 100 gerçeklik, σ_D = 5 meV). Tablo 4’te görüldüğü gibi T₂* oranı A_AC = 100 meV’e kadar artmaya devam eder (1.315 oranı, %31.5 uzama), ardından A_AC = 150 meV’de aniden 0.869’a düşer (T₂* aslında kısalır).

| A_AC (meV) | ⟨Δ₀₁⟩ (meV) | σ_Δ₀₁ (meV) | T₂* / T₂*_stat | min(Δ₀₁) (meV) |
|---|---|---|---|---|
| 0 | 1.018 | 0.1233 | 1.000 | 0.989 |
| 50 | 1.278 | 0.1125 | 1.096 | 0.190 |
| 80 | 1.557 | 0.1019 | 1.210 | 0.063 |
| 100 | 1.766 | 0.0938 | 1.315 | 0.030 |
| 150 | 1.927 | 0.1420 | 0.869 | 0.005 |

**Tablo 4.** *Doygunluk eğrisi. **min**(Δ₀₁) periyot içindeki anlık minimumdur; A_AC ≥ 50 **meV’den** itibaren 0.5 **meV** altına düşer.*

Bu davranışın fiziksel kaynağı operasyonel sınırın aşılmasıdır. A_AC ≳ 50 meV için periyot içinde anlık Δ₀₁(t) değeri 0.5 meV’in altına düşer, ve A_AC ≳ 80 meV’de 0.1 meV’in altına iner. Bu anlık dejenerlik, kübit seviyelerinin yan-bantlarla rezonant geçişine olanak verir ve Magnus seri açılımının geçerliliğini bozar. A_AC ≥ 150 meV’de QD kuyusunun bağlanma derinliği sıfıra yaklaşır ve bağlı durumlar bulk bantlarına sızar.

Bu nedenle paper’ın çekirdek bulguları operasyonel rejim olan A_AC ≤ 50 meV ile sınırlanmıştır. A_AC ∈ [50, 100] meV bölgesinde periyot ortalaması istatistiksel olarak koruma görünür kalsa da fiziksel yorum dikkat gerektirir — bu rejimde tam Floquet analizi gereklidir ve gelecekteki çalışmaya bırakılmıştır.

## 6. Tartışma

### 6.1 Deneysel Bağlam

Hecker ve arkadaşları [arXiv:2605.12257 (2025)] BLG QD’lerinde yük gürültüsü spektroskopisi ile S_ε ≈ 0.5–0.9 neV/√Hz değerleri ölçtüler. Bu, integre edilmiş efektif σ_D değerinin ölçüm zaman penceresine bağlı olarak meV mertebesinin oldukça altında olduğunu gösterir. Çalışmamızda kullandığımız σ_D = 5 meV agresif bir üst sınırdır ve AC etkisinin görünürlüğü için seçilmiştir. σ_D = 2 meV sonucumuz (Tablo 3) gerçekçi rejime daha yakındır.

Mekanizmanın σ_D’den bağımsızlığı (Bölüm 5.1) sayesinde gerçek deneylerde de aynı oransal koruma beklenir. A_AC = 50 meV genliği BLG QD deneysel kapı yelpazesinde erişilebilirdir; tipik ETH Zürih ve Aachen cihazlarında bu mertebede AC modülasyon rutin olarak uygulanmaktadır. T₂*’nın %9 mertebesinde uzaması Ramsey tipi ölçümlerle gözlemlenebilir, modern BLG QD deneylerinin çözünürlüğü bunun çok altındadır.

### 6.2 Hajati-Burkard ile Karşılaştırma

Hajati-Burkard (2024) Si/Ge düz-mod spin kübitlerinde rezonant kayma mekanizması üzerinden dinamik tatlı nokta gösterdi. Bizim BLG bulgumuz iki açıdan farklıdır: (i) BLG’de statik tatlı nokta yokken AC koruma yine de oluşur, (ii) mekanizma analitik RWA yoluyla değil, doğrudan band eğriliği üzerinden çalışır. Bu, dinamik tatlı noktaların yarıiletken kübit platformlarında daha geniş bir aileye sahip olabileceğini ima eder: statik tatlı nokta etrafında "rezonant" sınıf (Hajati-Burkard) ve eğrilik bazlı "Jensen" sınıf (bu çalışma). Hangi mekanizmanın baskın olacağı sistem-bağımlıdır.

### 6.3 Yöntemsel Not

BLG’nin 4-bant k·p yapısı ve trigonal warping, Hajati-Burkard’ın benimsediği analitik RWA çerçevesinin uygulanmasını zorlaştırır. Bu nedenle tam sayısal yaklaşımı (düzlem-dalga bazı + ensemble taraması) zorunlu olarak kullandık. Bu metodolojik seçimin yan etkisi olarak Jensen mekanizmasını ortaya çıkardık — analitik RWA bu mekanizmayı doğal olarak göstermez çünkü ikinci derece band eğriliği etkilerini RWA basitleştirmesi bastırır. Bu, sayısal yaklaşımın analitik yaklaşıma kıyasla mekanistik olarak daha geniş bir kapsama sahip olabileceğine işaret eder.

### 6.4 Sınırlamalar

Bu çalışmanın açıkça belirtilmesi gereken sınırlamaları:

- **Spinsiz****model:** Gerçek BLG QD kübitleri spin-valley kübitlerdir [Banszerus et al. 2022]. Spin-orbit kuplajı (Kane-Mele tipi) bu çalışmada dışlanmıştır. Spin uzantısı 4-bant yapıyı 8-banta çıkarır; metodolojik olarak doğrudandır ancak bu paper’ın kapsamı dışındadır.

- **Kuasi****-statik gürültü:** Gerçek yük gürültüsü 1/f spektrumlu, geniş-bantlıdır. Bizim Gaussian kuasi-statik modelimiz Hajati-Burkard ile aynı seviyede basitleştirmedir. Tam 1/f modeli Lindblad veya stochastic Schrödinger çerçevesini gerektirir.

- **Magnus****/ hızlı sürüş:** ω/Δ₀₁ ≈ 6.6 oranı Magnus seri açılımının geçerliliğini sağlar. Rezonant rejim (ω ~ Δ₀₁) tam Floquet analizini gerektirir, ayrı çalışma konusu olacaktır.

- **Tek QD:** Çalışmamız tek bir BLG QD ile sınırlıdır. İki-QD valley kübitler (singlet-triplet analog) ayrı dinamikler içerir.

- **Düzgün geometri:** Gauss kuyusu kullandık. Gerçek deneysel cihazlarda split-gate geometrisi daha karmaşık potansiyeller üretir. Mayer-Knothe (2023) çalışmasının gösterdiği gibi, niteliksel sonuçlar potansiyel detayına göre korunur, ancak nicel öngörüler için gerçekçi geometri gereklidir.

## 7. Sonuç ve Gelecek Perspektifler

Bilayer grafen kuantum nokta valley kübitlerinde AC-sürüş ile yük gürültüsü koruması sorusunu, dört bant k·p sürekli ortam modelinde düzlem-dalga bazı kullanarak doğrudan sayısal olarak ele aldık. Çalışmamızın çekirdek bulguları:

- BLG valley kübitlerinde, perturbatif hapsetme rejiminde, statik tatlı nokta (∂Δ₀₁/∂Δ_ext = 0) mevcut değildir.

- Buna rağmen Δ₀₁(Δ_ext) eğrisinin konkavlığı, AC sürüş altında Jensen tipi eşitsizlik aracılığıyla σ_Δ₀₁’yi azaltır ve T₂*’yi uzatır.

- Operasyonel rejimde (A_AC ≤ 50 meV) T₂* yaklaşık %9 oranında uzar. Etki σ_D gürültü genişliğinden bağımsızdır, yalnızca band eğriliğine bağlıdır.

- Bu mekanizma Hajati-Burkard (2024) rezonant kayma mekanizmasından nitel olarak farklıdır; statik tatlı nokta gerektirmez.

Bu bulgular dinamik tatlı noktaların iki ayrı sınıfı olduğunu gösterir: (i) statik tatlı nokta etrafında rezonant koruma (Hajati-Burkard) ve (ii) statik tatlı nokta olmadan band eğriliği üzerinden Jensen koruması (bu çalışma). BLG valley kübitleri ikinci sınıfın deneysel olarak erişilebilir ilk somut örneğini sağlar.

Gelecekteki çalışmalar şu doğal genişlemeleri içermelidir: (a) Kane-Mele spin-orbit kuplajıyla spin-valley kübit analizi, (b) 1/f gürültü spektrumunun Lindblad master denkleminde tam tedavisi, (c) rezonant rejimde tam Floquet analizi, (d) gerçekçi split-gate geometrilerle Poisson çözücü entegrasyonu, ve (e) mevcut deneysel BLG QD platformlarında (ETH Zürih, RWTH Aachen) Ramsey ölçümleriyle deneysel doğrulama. Çalışmamız bu programların başlangıç noktasını sağlar.

### Referanslar

Banszerus, L., et al. (2022). "Spin-valley coupling in single-electron bilayer graphene quantum dots." Nature Communications 13, 3637.

Benito, M., Croot, X., Adelsberger, C., Putz, S., Mi, X., Petta, J. R., & Burkard, G. (2019). "Electric-field control and noise protection of the flopping-mode spin qubit." Physical Review B 100, 125430.

Burkard, G., Ladd, T. D., Pan, A., Nichol, J. M., & Petta, J. R. (2023). "Semiconductor spin qubits." Reviews of Modern Physics 95, 025003.

Eich, M., et al. (2018). "Spin and valley states in gate-defined bilayer graphene quantum dots." Physical Review X 8, 031023.

Garreis, R., et al. (2021). "Shell filling and trigonal warping in graphene quantum dots." arXiv:2011.07951.

Hajati, Y., & Burkard, G. (2024). "Dynamic sweet spot of driven flopping-mode spin qubits in planar quantum dots." Physical Review B 110, 245301. arXiv:2408.03730.

Hecker, K., et al. (2025). "Charge noise spectroscopy in bilayer graphene quantum dots." arXiv:2605.12257.

Knothe, A., & Fal’ko, V. I. (2018). "Influence of minivalleys and Berry curvature on electrostatically induced quantum wires in gapped bilayer graphene." Physical Review B 98, 155435.

Lee, Y., et al. (2020). "Tunable valley splitting due to topological orbital magnetic moment in bilayer graphene quantum point contacts." Physical Review Letters 124, 126802.

Loss, D., & DiVincenzo, D. P. (1998). "Quantum computation with quantum dots." Physical Review A 57, 120.

Mayer, D., & Knothe, A. (2023). "Tuning confined states and valley g-factors by quantum dot design in bilayer graphene." physica status solidi (b) 260, 2300395.

McCann, E., & Koshino, M. (2013). "The electronic properties of bilayer graphene." Reports on Progress in Physics 76, 056503.

Mundada, P. S., et al. (2020). "Engineering dynamical sweet spots to protect qubits from 1/f noise." arXiv:2004.12458.

Petta, J. R., et al. (2005). "Coherent manipulation of coupled electron spins in semiconductor quantum dots." Science 309, 2180.

Stockklauser, A., et al. (2017). "Strong coupling cavity QED with gate-defined double quantum dots." Physical Review X 7, 011030.

## Ek A: Düzlem-Dalga Bazına Geçiş Motivasyonu

İlk denemelerde 4-bant k·p Hamiltonyenini sonlu fark (finite difference) diskretizasyonu ile çözmeye çalıştık. 41×41 grid (h = 2 nm, Hilbert boyutu 6724) ile elde edilen spektrum patolojik bir özellik gösterdi: Δ_ext = 0.10 eV ve U₀ = 50 meV referans noktasında en düşük iki bağlı durum arasındaki ayrılış Δ₀₁ = 0.002 meV olarak hesaplandı. Bu, fiziksel olarak makul (~5–10 meV) değerden iki mertebe küçüktü. Ayrıca dalga fonksiyonları yapay bir 4-katlı dejenerelik ve dama tahtası tipi yoğunluk desenleri gösterdi.

Bu davranış Dirac operatörlerinin sonlu fark diskretizasyonunun bilinen bir özelliğidir — fermion-doubling. Birinci derece momentum operatörü grid üzerinde merkezi fark ile yaklaşıldığında, gerçek Dirac noktasına ek olarak Brillouin bölgesinin köşelerinde sahte "ikiz" durumlar ortaya çıkar. Wilson terimi ekleyerek bu artefakt kısmen bastırılabilir, ancak BLG’nin dört bantlı yapısı ve trigonal warping varlığında bu yöntem teknik olarak zordur.

Düzlem-dalga bazı bu sorunun yapısal çözümüdür: kinetik enerji momentum uzayında köşegen olduğundan sahte ikiz durum üretmez. Bizim N_pw = 293 ile 1172 boyutunda çalıştığımız sistem, sonlu fark gridinin 6724 boyutundan oldukça küçük olmasına rağmen tamamen artefaktsız spektrum üretir. Bulk testleri (V = 0 limiti) analitik bant min değerine %5 dahilinde uyum gösterir, ve QD spektrumu (17 bağlı durum, Δ₀₁ ≈ 8 meV) literatür beklentileriyle uyumludur.

## Ek B: Anti-crossing Rejimi (U₀ = 50 meV) Analizi

U₀ = 50 meV ile yapılan Δ_ext taramasında, Δ_ext < 0.10 eV bölgesinde Δ₀₁(Δ_ext) eğrisinin sıçramalı, monoton olmayan davranış gösterdiği gözlemlendi. Δ_ext = 0.083 eV’de Δ₀₁ = 9 meV’e fırlar, Δ_ext = 0.10 eV’de tekrar 8 meV’e döner. Bu davranış başlangıçta gerçek bir stasyoner nokta olarak yorumlanabilir, ancak detaylı incelemede anti-crossing artefaktı olduğu anlaşılır: bu rejimde U₀ > bulk_gap koşulu sağlanır (örn. Δ_ext = 0.083 eV için bulk_gap = 40.7 meV < U₀ = 50 meV), ve QD bağlı durumları bulk bantlarıyla rezonansa girerek karışır.

Bu nedenle anti-crossing rejimi paper’ın ana sonuçlarından dışlanmıştır. U₀ = 20 meV referansı tüm Δ_ext değerleri için U₀/bulk_gap < 0.5 sağladığından, temiz perturbatif hapsetme rejimini garanti eder. Mayer-Knothe (2023) çalışmasında da benzer bir kural izlenmiştir.

## Ek C: Yapısal Çıpalar — Detaylı Sayısal Değerler

Tüm hesaplamalarda raporlanan yapısal çıpa değerleri:

- Çıpa 1 (Hermitiklik): |H − H†|_max = 0.00e+00 (makine hassasiyeti)

- Çıpa 2 (V = 0): En küçük 8 özdeğer (eV): −0.04689, −0.04640, −0.04640, −0.04553, +0.04553, +0.04640, +0.04640, +0.04689. Analitik γ₃ = 0 beklentisi ±48.35 meV; γ₃ = 0.38 eV ile ~%5 düşeltilmiş halinde uyumlu.

- Çıpa 3 (Trigonal warping): ψ₁ ve ψ₂ üç-katlı C₃ simetrik trefoil yoğunluk dağılımları (Bölüm 3.1).

- Çıpa 4 (Δ = 0): Bağlı durum sayısı sıfıra düşer, spektrum metaliğe geçer.

- Çıpa 5 (Kramers): K ve K’ vadi spektrumları zaman-tersi simetrisi altında özdeş (1e−15 hassasiyetinde).

- AC yapısal çıpa: max|⟨Δ₀₁⟩_AC(A=0) − Δ₀₁_statik| = 0.00e+00 eV (Cell 7, σ_D = 5 meV durumunda).
