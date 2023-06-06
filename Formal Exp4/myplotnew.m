semilogy(1:100:Itrmax1,err_SASC_hist(1:100:Itrmax1),1:100:Itrmax1,err_nested_SGD_hist(1:100:Itrmax1)...
     );
legend('SASC-SGD','Nested-SGD');
xlabel('Number of Stochastic Evaluations');
ylabel('Relative Error');
ylim([1e-4 10]);