# BÁO CÁO PHÂN TÍCH VÀ TRÌNH BÀY DỰ ÁN PRISM (PLOTS/PRISM)
## Privacy-preserving Rare Disease Analysis using Fully Homomorphic Encryption

---

### TỔNG QUAN VỀ DỰ ÁN PRISM

PRISM là một giải pháp mật mã tiên tiến giải quyết bài toán **phân tích biến thể gen tìm bệnh hiếm bảo vệ quyền riêng tư (Privacy-preserving Rare Disease Analysis)** sử dụng **Mã hóa hoàn toàn đồng hình (Fully Homomorphic Encryption - FHE)** và **Mật mã ngưỡng nhiều bên (Multi-party Threshold FHE)**.

Trong y khoa và di truyền học, việc phát hiện các biến thể gen gây ra bệnh hiếm đòi hỏi phải so sánh mẫu ADN của người bệnh (cases) và người khỏe mạnh (controls) hoặc giữa các thành viên trong gia đình (bố, mẹ, con). Tuy nhiên, dữ liệu gen chứa thông tin cực kỳ nhạy cảm về nhận dạng và y sinh. PRISM cho phép các trung tâm y tế hợp tác phân tích dữ liệu gen trên quy mô lớn mà **không cần chia sẻ dữ liệu thô (raw genomic data)** và **không lộ thông tin cá nhân của bệnh nhân**.

---

### 1. TÍNH MỚI VÀ ĐÓNG GÓP NỔI BẬT (NOVELTY)

PRISM giải quyết các rào cản lớn nhất của FHE trong phân tích gen (độ trễ cao, giới hạn chiều sâu phép nhân, chi phí giải mã cao):

1. **Thuật toán lọc biến thể tối ưu chiều sâu phép nhân (Add-in Method):**
   - Các phương pháp FHE truyền thống dùng phép nhân đồng hình $O(M)$ để kiểm tra điều kiện lọc trên $M$ mẫu, dẫn đến **Multiplicative Depth** rất cao ($M \ge 12$). Điều này đòi hỏi tham số FHE rất lớn (kích thước ciphertexts khổng lồ, tính toán rất chậm).
   - PRISM đề xuất thuật toán **Add-in**: thay thế phép nhân bằng phép cộng đồng hình chuỗi bit biến thể và phép trừ ngưỡng cố định, giảm chiều sâu phép nhân xuống **độ sâu cố định = 2** (hoặc = 1) độc lập với số lượng mẫu $M$.

2. **Cơ chế Mã hóa đồng hình Threshold Multi-party (Multiparty FHE - MFHE):**
   - Không phụ thuộc vào một bên thứ ba tin cậy (Trusted Third Party - TTP).
   - Dữ liệu được mã hóa bằng **khóa công khai chung (Collective Public Key)** do $P$ bên cùng tạo ra.
   - Quá trình giải mã kết quả lọc chỉ thành công khi có sự đồng thuận của tất cả $P$ bên tham gia (sử dụng Partial Decryption & Fusion).

3. **Cơ chế Làm mù ngẫu nhiên (Random Blinding) bảo mật tuyệt đối kết quả trung gian:**
   - Để ngăn chặn việc bên tính toán (Cloud server) hoặc các bên tham gia rút trích thông tin từ vị trí biến thể trùng khớp, PRISM nhân vector kết quả mã hóa với một vector số ngẫu nhiên bí mật trước khi xáo trộn (shuffle) và giải mã.
   - Nhờ đó, server giải mã chỉ biết **tổng số lượng biến thể thỏa mãn** mà **không biết vị trí cụ thể** của các biến thể gen đó.

4. **Tối ưu hóa mã hóa vector SIMD (Batching / Packing):**
   - Sử dụng kỹ thuật BFV SIMD Packing, đóng gói hàng chục ngàn biến thể gen ($N = 25,000$) vào trong **một Ciphertext duy nhất**, giúp giảm thời gian tính toán và bộ nhớ lưu trữ đi hàng nghìn lần.

---

### 2. MÔ HÌNH DI TRUYỀN VÀ THUẬT TOÁN LỌC (INHERITANCE MODELS & ALGORITHMS)

PRISM hỗ trợ 2 mô hình di truyền bệnh hiếm phổ biến nhất:

#### A. Mô hình Di truyền Lặn & Trội (Recessive and Dominant Inheritance Models)
Mô hình này so sánh cấu trúc kiểu gen của mẫu với một truy vấn định sẵn (Query). Mỗi vị trí biến thể gen của một mẫu được mã hóa bởi 2 bit $(b_0, b_1)$ đại diện cho 2 alen (00: Homozygous Reference, 01/10: Heterozygous, 11: Homozygous Alternate).

- **Phép toán kiểm tra khớp bit (Encrypted Bit Matching - XOR Gate):**
  Cho bit dữ liệu $c$ và bit truy vấn $q$:
  $$\text{Match}(c, q) = 1 - (c - q)^2 = 1 - (c - q)(c - q)$$
  Nếu $c = q$, $\text{Match}(c, q) = 1$. Ngược lại bằng 0.

- **Thuật toán Mul-in (Phương pháp Nhân):**
  Để biến thể $i$ thỏa mãn truy vấn trên tất cả $2M$ bit của $M$ mẫu:
  $$\text{Result}_i = \prod_{j=1}^{2M} \text{Match}(c_{j, i}, q_j)$$
  *Chiều sâu phép nhân:* $O(\log_2(2M))$, yêu cầu $MultDepth \ge 12$ cho 16 mẫu.

- **Thuật toán Add-in (Phương pháp Cộng - Đóng góp cốt lõi của PRISM):**
  PRISM cộng tất cả các kết quả match lại:
  $$\text{Sum}_i = \sum_{j=1}^{2M} \text{Match}(c_{j, i}, q_j)$$
  Nếu biến thể $i$ khớp toàn bộ $2M$ bit, thì $\text{Sum}_i = 2M$.
  PRISM thực hiện phép trừ ngưỡng:
  $$\text{Diff}_i = \text{Sum}_i - 2M$$
  Nếu biến thể khớp hoàn toàn, $\text{Diff}_i = 0$. Nếu không khớp, $\text{Diff}_i < 0$ (khác 0).

  Sau đó nhân với vector ngẫu nhiên $R_i \neq 0$:
  $$\text{Blind}_i = \text{Diff}_i \times R_i$$
  - Nếu khớp hoàn toàn: $0 \times R_i = 0$.
  - Nếu không khớp: $\text{Diff}_i \times R_i \neq 0$.
  *Chiều sâu phép nhân:* Chỉ bằng **2**!

---

#### B. Mô hình Di truyền Mới phát sinh (De Novo Mutation Model)
Mô hình De Novo tìm kiếm đột biến xuất hiện ở con (Proband) nhưng không có ở Bố và Mẹ (Parents).

- Quy tắc: Mẹ ($M$), Bố ($F$), Con ($C$).
  - Biến thể De Novo xảy ra khi Con mang alen đột biến ($C_j = 1$) nhưng Mẹ và Bố không mang alen đó ($M_j = 0$ và $F_j = 0$).

- **Thuật toán Add-in cho De Novo:**
  Đối với $M$ mẫu gia đình (Trio):
  PRISM tính toán:
  $$\text{Sum}_i^{(0)} = \sum_{\text{Con}} C_i + \sum_{\text{Bố/Mẹ}} (1 - P_i)$$
  Nếu đúng chuẩn De Novo tại biến thể $i$, $\text{Sum}_i^{(0)} = M$.
  Biến đổi: $\text{Diff}_i = \text{Sum}_i^{(0)} - M$.
  Nhân với số ngẫu nhiên $R_i$: nếu $\text{Diff}_i = 0$ thì kết quả bằng 0 (khớp).

---

### 3. KIẾN TRÚC VÀ CÁCH HOẠT ĐỘNG TOÀN DIỆN (WORKFLOW)

1. **Giai đoạn Key Generation (Thiết lập khóa Multi-party):**
   - $P$ phòng thí nghiệm / bệnh viện phối hợp tạo ra Khóa công khai chung $PK_{\text{joint}}$ và các khóa hỗ trợ tính toán đồng hình ($KeyEvalMult$, $KeyEvalSum$).
   - Mỗi bên giữ một phần khóa bí mật $SK_p$.

2. **Giai đoạn Mã hóa & Đóng gói dữ liệu (Encoding & Encryption):**
   - Dữ liệu vcf/genotype được chia thành từng khối $N = 25,000$ biến thể.
   - Mỗi khối được mã hóa bằng BFV Scheme với $PK_{\text{joint}}$ thành Ciphertext $CT$.

3. **Giai đoạn Tính toán đồng hình trên Cloud (Encrypted Filtering):**
   - Cloud Server nhận các Ciphertext từ các bên.
   - Chạy thuật toán Lọc biến thể (Add-in hoặc Mul-in) hoàn toàn trên miền mã hóa.
   - Áp dụng Random Blinding và Shuffle vị trí các Ciphertext.

4. **Giai đoạn Giải mã ngưỡng (Threshold Decryption):**
   - Kết quả mã hóa được gửi tới $P$ bên.
   - Mỗi bên dùng $SK_p$ thực hiện giải mã một phần (Partial Decryption).
   - Tổng hợp các bản giải mã một phần (Fusion) để thu được kết quả cuối cùng: số lượng biến thể gen khớp với tiêu chuẩn chẩn đoán.

---

### 4. BẢNG SO SÁNH ADD-IN VS MUL-IN

| Tiêu chí | Phương pháp Mul-in (Cũ) | Phương pháp Add-in (PRISM) |
| :--- | :--- | :--- |
| **Chiều sâu phép nhân ($MultDepth$)** | Rất cao ($12 - 16$) | Cực thấp ($1 - 2$) |
| **Kích thước Ciphertext** | Cực lớn ($\sim 10 - 50$ MB/CT) | Nhỏ ($\sim 1 - 2$ MB/CT) |
| **Thời gian lọc trên Cloud** | Chậm (hàng chục phút) | Cực nhanh (vài giây đến vài phút) |
| **Bảo mật vị trí biến thể** | Phụ thuộc vào giải mã tổng | Dùng Random Blinding + Shuffle |
| **Khả năng mở rộng số mẫu ($M$)** | Khó mở rộng ($MultDepth$ tăng lên) | Mở rộng tuyệt vời ($MultDepth$ luôn = 2) |

---

### 5. CHUYỂN ĐỔI SANG JAVASCRIPT & XÂY DỰNG DEMO

Để minh họa trực quan thuật toán cho GS, chúng ta xây dựng hệ thống bằng **JavaScript / Node.js**:
1. **Core FHE Simulation & Engine (`prism-engine.js`):** Mô phỏng chính xác thuật toán mã hóa BFV, mã hóa vector SIMD, các cổng XOR đồng hình, thuật toán Add-in & Mul-in, cơ chế Blinding và Giải mã Multi-party.
2. **CLI Benchmark Script (`prism-demo.js`):** Chạy thử nghiệm tự động trên tập dữ liệu mẫu, so sánh chính xác kết quả giữa tính toán rõ (Cleartext Ground Truth) và tính toán FHE mã hóa.
3. **Interactive Web UI Dashboard (`web-demo/index.html`):** Giao diện tương tác trực quan hiển thị dữ liệu gen, từng bước tính toán FHE, cho phép thay đổi tham số $M$ (số mẫu), $N$ (số biến thể) và kiểm tra kết quả tức thì.
