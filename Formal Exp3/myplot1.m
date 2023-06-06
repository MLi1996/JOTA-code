semilogy(1:100:Itrmax1,err_SASC_hist(1:100:Itrmax1),1:100:Itrmax1,err_nested_SGD_hist(1:100:Itrmax1),...
     1:100:Itrmax,err_nested_SGDM_hist(1:100:Itrmax),1:100:Itrmax,err_static_SGDM_hist(1:100:Itrmax));
legend('SASC-SGD','Nested-SGD','Nested-SGDM','Static-SGDM');
xlabel('Number of Stochastic Evaluations');
ylabel('Relative Error');
ylim([1e-4 10]);