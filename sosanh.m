% 1. Đọc dữ liệu từ tệp CSV
% MATLAB sẽ tự động xử lý khoảng trắng trong tên cột (ví dụ: 'centrealize 2' thành 'centrealize2')
% Thay đường dẫn này bằng đường dẫn thực tế trên máy bạn
fullpath = 'D:\báo\bait2\cod\compare.xlsx'; 

opts = detectImportOptions(fullpath);
data = readtable(fullpath, opts);

% 2. Hiển thị thông tin bảng để kiểm tra tên cột
disp('Thông tin dữ liệu:');
disp(head(data));

% 3. Vẽ biểu đồ
figure;
hold on; % Giữ đồ thị để vẽ nhiều đường trên cùng một trục

% Vẽ cột 'distributed'
plot(data{:, 1}, '-o', 'LineWidth', 1.5, 'DisplayName', 'Distributed');

% Vẽ cột 'centrealize 2' (Lưu ý: readtable thường đổi tên có khoảng trắng thành 'centrealize2')
% Nếu gặp lỗi, hãy kiểm tra lại tên cột bằng lệnh data.Properties.VariableNames
plot(data{:, 2}, '-s', 'LineWidth', 1.5, 'DisplayName', 'Centralized');

% 4. Trang trí biểu đồ
grid on;
xlabel('Iter');
ylabel('EE (Mbits/J)');
title('So sánh Distributed và Centralized');
legend('show', 'Location', 'best');

% Tùy chỉnh thêm (nếu cần)
set(gca, 'FontSize', 12);

hold off;