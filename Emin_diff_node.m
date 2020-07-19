filename = 'A_structural_connectivity_matrix.xlsx';
A = xlsread(filename);

filename = 'B_control_input_matrix.xlsx';
B = xlsread(filename);

x0 = [8.225387253	5.563230023	7.403526234	4.396248436	12.43532793	12.34020411	7.99205084	4.052766589	7.348044999	12.48262293	6.529798392	1.987671184	1.337519158	2.443223354	2.176707023	2.963241306	2.930496026	1.832843636	3.180074586];
xf = [26.60157269	13.30354717	5.840641061	3.755314359	12.52854331	14.83160509	10.24569555	3.772343079	10.2748911	11.14400986	9.526875948	2.744283669	2.304353076	6.429682189	2.951137898	5.880818184	8.071683663	4.200881174	14.70349171];

T = 1

nor = 1;


[ x, u, n_err ] = min_eng_cont(A, T, B, x0', xf', nor)

Emin_all = sum((sum(u.^2,1).*(T/1000)))
%loop function to calculate increase in min control energy with each node
%suppressed
     for node = 1:19
     B = eye(19)
     B(node, node) = 0
     [ x, u, n_err ] = min_eng_cont(A, T, B, x0', xf', nor)
     Emin(node) = sum((sum(u.^2,1).*(T/1000)))
     Emin_dif(node) = Emin(node) - Emin_all
end