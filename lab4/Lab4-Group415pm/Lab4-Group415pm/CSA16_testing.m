%% Before running, set up the testbench cell name, input and output table, clock period etc. In the "EDIT BELOW/ABOVE ACCORDING TO YOUR DUT" block
% (c) 2019 Oscar Castaneda, Olalekan Afuye, Charles Jeon & Christoph Studer
% Modified by Yixiao Du (yd383@cornell.edu)
% ECE4740 logic module test tool.

clear


% vvvvvvvvvvvvvvvvvvvvvvvvvvv EDIT BELOW ACCORDING TO YOUR DUT vvvvvvvvvvvvvvvvvvvvvvvvvvv
% set up the name and parameters of your testbench cell
tb_name = 'lab4_tb';
% set up inputs and outputs, MSB is on the left, use ... to contiune a line
input_nbr = 33;
input_table = {
  'A\<15\>','A\<14\>','A\<13\>','A\<12\>','A\<11\>','A\<10\>','A\<9\>','A\<8\>','A\<7\>','A\<6\>','A\<5\>','A\<4\>','A\<3\>','A\<2\>','A\<1\>','A\<0\>',...
  'B\<15\>','B\<14\>','B\<13\>','B\<12\>','B\<11\>','B\<10\>','B\<9\>','B\<8\>','B\<7\>','B\<6\>','B\<5\>','B\<4\>','B\<3\>','B\<2\>','B\<1\>','B\<0\>',...
  'Cin'
};
output_nbr = 17;
output_table = {
  'Cout',...
  'S\<15\>','S\<14\>','S\<13\>','S\<12\>','S\<11\>','S\<10\>','S\<9\>','S\<8\>','S\<7\>','S\<6\>','S\<5\>','S\<4\>','S\<3\>','S\<2\>','S\<1\>','S\<0\>',
};

% define total test case numbers
sample_nbr = 10; 
% define period (in ps)
cp  =  2000;
% define rise and fall time (in ps)
t_rf = 10;
% initial setup time
t_ini = cp;
% maximum propagation delay (in ps)
max_tp = 1300;
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^ EDIT ABOVE ACCORDING TO YOUR DUT ^^^^^^^^^^^^^^^^^^^^^^^^^^^


% set up cds_srr function
addpath('/opt/cadence/MMSIM151/tools.lnx86/spectre/matlab/64bit');

% directory that contains the simulation outputs
directory = sprintf('%s/Cadence/simulation/%s/spectre/schematic/psf', getenv('HOME'), tb_name);

% set up basic parameters
Vdd = 1.2;

% total simulation data length
a = cds_srr(directory, 'tran-tran', input_table{1}, 0);
data_len = max(size(a.V));

% get input signals
input_mtx = zeros(data_len,input_nbr);
for i = 1:input_nbr
  a = cds_srr(directory, 'tran-tran', input_table{i}, 0); 
  % extract voltages of signals
  input_mtx(:,i) = a.V;
end

% convert time into ps
t_ps = a.time*1e12;

% get output signals
output_mtx = [];
for i=1:output_nbr
  y = cds_srr(directory, 'tran-tran', output_table{i}, 0);
  output_mtx = [output_mtx y.V];
end

% we sample the inputs at the middle of a cycle
t_ps_sample_in = t_ini + cp/2 + (0:sample_nbr-1)*cp;

% we sample the outputs tp_max after an input changes (each cp),
t_ps_sample_out = t_ini + max_tp + (0:sample_nbr-1)*cp;

%% DUT output
% create base for expected output waveform
input_vec_bits = [];
for i = 1:input_nbr
  input_vec_bits = [input_vec_bits (input_mtx(:,i) > (Vdd/2))];  
end


% vvvvvvvvvvvvvvvvvvvvvvvvvvv EDIT BELOW ACCORDING TO YOUR DUT vvvvvvvvvvvvvvvvvvvvvvvvvvv
% define the functionality of DUT
input_A_dec = bi2de(input_vec_bits(:,(1:16)),'left-msb');
input_B_dec = bi2de(input_vec_bits(:,(17:32)),'left-msb');
input_Cin_dec = bi2de(input_vec_bits(:,33),'left-msb');
exp_dec = input_A_dec + input_B_dec + input_Cin_dec;
s_expected  = mod(exp_dec,(2^16));
cout_expected = floor(exp_dec./(2^16));
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^ EDIT ABOVE ACCORDING TO YOUR DUT ^^^^^^^^^^^^^^^^^^^^^^^^^^^


%Check each one of the sampling points
err_flag = 0;
dut_output = zeros(sample_nbr,output_nbr);
sample_pt = t_ps_sample_out;
disp('A + B + Cin = Cout,S')
for i=1:sample_nbr
    % find t_ps closest (from the right) to the t_ps_sample_in and _out
    t_ps_idx_in  = find(t_ps-t_ps_sample_in(i)>=0,1);
    t_ps_idx_out = find(t_ps-t_ps_sample_out(i)>=0,1);
    
    % measure the outputs and declare 1 if it is greater than Vdd/2    
    dut_output(i,:) = output_mtx(t_ps_idx_out,:) > (Vdd/2);


    % vvvvvvvvvvvvvvvvvvvvvvvvvvv EDIT BELOW ACCORDING TO YOUR DUT vvvvvvvvvvvvvvvvvvvvvvvvvvv
    % interpret UDT outputs
    s_result  = bi2de(dut_output(:,2:17),'left-msb');
    cout_result = bi2de(dut_output(:,1),'left-msb');
    
    % check results
    err = (cout_expected(t_ps_idx_out,:) ~= cout_result(i,:)) || (s_expected(t_ps_idx_out,:) ~= s_result(i,:));
    if (err)
      disp('---------------------------------------------------')
      disp('[ERROR!] Expected:')
      msg = sprintf('%x + %x + %x = %x,%x',input_A_dec(t_ps_idx_in),input_B_dec(t_ps_idx_in),...
                    input_Cin_dec(t_ps_idx_in),cout_expected(t_ps_idx_out),s_expected(t_ps_idx_out));
      disp(msg)
      disp('But actually get:')
      msg = sprintf('%x + %x + %x = %x,%x',input_A_dec(t_ps_idx_in),input_B_dec(t_ps_idx_in),...
                    input_Cin_dec(t_ps_idx_in),cout_result(i),s_result(i));
      disp(msg)
      err_flag  = err_flag + 1;
    else
      disp('---------------------------------------------------')
      msg = sprintf('%x + %x + %x = %x,%x',input_A_dec(t_ps_idx_in),input_B_dec(t_ps_idx_in),...
                    input_Cin_dec(t_ps_idx_in),cout_result(i),s_result(i));
      disp(msg)
    end  
    % ^^^^^^^^^^^^^^^^^^^^^^^^^^^ EDIT ABOVE ACCORDING TO YOUR DUT ^^^^^^^^^^^^^^^^^^^^^^^^^^^


end
disp('---------------------------------------------------')
if err_flag == 0
    disp('Wryyyy! Your circuit has no errors!')
    disp('---------------------------------------------------')
end