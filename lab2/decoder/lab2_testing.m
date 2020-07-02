%% Before running, set up the testbench cell name on line 8
% your decoder outputs should be labeled as Y0,Y1,Y2,Y3
% your decoder inputs should be labeled as A,B
%
% (c) 2019 Oscar Castaneda, Olalekan Afuye, Charles Jeon & Christoph Studer

% set up the name of your testbench cell
tb_name = 'lab2_tb';

% set up cds_srr function
addpath('/opt/cadence/MMSIM151/tools.lnx86/spectre/matlab/64bit');

% directory that contains the simulation outputs
directory = sprintf('%s/Cadence/simulation/%s/spectre/schematic/psf', getenv('HOME'), tb_name);

% set up basic parameters
Vdd = 1.2; % define vdd

% define period (in ps)
period_a = 1000; % A
period_b =  500; % B

% get input signals
a = cds_srr(directory, 'tran-tran', 'A', 0);
b = cds_srr(directory, 'tran-tran', 'B', 0);

% convert time into ps
t_ps = a.time*1e12;

% extract voltages of signals
a = a.V;
b = b.V;

% get output signals and put them together in a table where the i-th
% column corresponds to the 'Y(i-1)' output
y_mtx = [];
for i=1:4
    signal_name = ['Y',int2str(i-1)];
    y = cds_srr(directory, 'tran-tran', signal_name, 0);
    y_mtx = [y_mtx y.V];
end

exp_y_mtx = zeros(size(y_mtx));
sample_wvf = zeros(size(y_mtx));
mydecoder_output = zeros(4,4);
exp_decoder_output = zeros(4,4);

% we sample the inputs at the middle of a cycle
t_ps_sample_in = 6*period_a + period_b/2 + (0:3)*period_b;

% we sample the outputs 230ps after an input changes (each 500ps),
% during the fourth time the inputs repeat
t_ps_sample_out = 6*period_a + 10 + 230 + (0:3)*period_b;

%% decoder output

% create base for expected output waveform
a_bits = (a > Vdd/2);
b_bits = (b > Vdd/2);
vec_bits = [a_bits b_bits];
exp_dec = bi2de(vec_bits,'left-msb');

%Check each one of the sampling points
err_flag = 0;
for i=1:4
    % find t_ps closest (from the right) to the t_ps_sample_in and _out
    t_ps_idx_in  = find(t_ps-t_ps_sample_in(i)>=0,1);
    t_ps_idx_out = find(t_ps-t_ps_sample_out(i)>=0,1);
    
    % measure the outputs and declare 1 if it is greater than Vdd/2    
    mydecoder_output(i,:) = y_mtx(t_ps_idx_out,:) > (Vdd/2);
    
    %create expected output waveform
    exp_y_mtx(:,i) = Vdd*(exp_dec == (i-1));
    
    %create sampling waveform
    sample_wvf(t_ps_idx_out,exp_dec(t_ps_idx_in)+1) = Vdd;
    
    % expected decoder output is given by:
    exp_decoder_output(i,exp_dec(t_ps_idx_in)+1) = 1;
    
    if (sum(exp_decoder_output(i,:) ~= mydecoder_output(i,:)))

    disp(['Expected output for input '...
        'A=' num2str(vec_bits(t_ps_idx_in,1)) ...
        ' B=' num2str(vec_bits(t_ps_idx_in,2)) ...        
        ' is y' num2str(4-i) ...
        ' but measured output is y' num2str(exp_dec(t_ps_idx_out))...
        ]) 
    err_flag  = err_flag + 1;
    end
end

if err_flag == 0
    disp('Your decoder circuit has no errors :)')
end
    %
%% plots
figure(1)
set(gcf,'units','pixel');
set(gcf,'position',[0,200,800,600]);
% we have 4 outputs so total time period elapsed is 2*period_a = 2000ps
subplot(2,1,1);
plot(t_ps,a,'b',t_ps,b,'r:','linewidth',2)
legend('A','B')
grid on
% note that our simulations start after 5ns so we add it to our window
xlim(6*period_a + [0 2*period_a])
title('inputs')
xlabel('time [ps]')
ylabel('input to decoder [V]')
set(gca,'fontsize',14)

% output from your circuit
subplot(2,1,2);
plot(t_ps,y_mtx,'linewidth',2)
legend('Y0','Y1','Y2','Y3')
grid on
xlim(6*period_a + [0 2*period_a])
title('outputs')
xlabel('time [ps]')
ylabel('output [V]')
set(gca,'fontsize',14)

% print to SVG (vector graphics file)
print -dsvg Lab2_decoder_IO.svg

figure(2)
set(gcf,'units','pixel');
set(gcf,'position',[800,200,800,600]);
for i=1:4
    subplot(2,2,i)
    plot(t_ps,y_mtx(:,i),'k',t_ps,exp_y_mtx(:,i),'b:',t_ps,sample_wvf(:,i),'r--','linewidth',3);
    grid on
    legend('actual output','ideal output','sampling point','location','south')
    xlim(6*period_a + [0 2*period_a])
    ylim([-1 1.5])
    xlabel('time [ps]')
    ylabel('output [V]')
    title(['y' num2str(i-1)])
    set(gca,'fontsize',10)
end

% print to SVG (vector graphics file)
print -dsvg Lab2_decoder_expresp.svg