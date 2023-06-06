Itr2_SVRG = ceil(Itrmax2*3/4);
Itr2_SS = ceil(Itrmax2*3/4);
Itr2_Catalyst = ceil(Itrmax2^2/ritr_num_Catalyst);
Itr2_CS = ceil(Itrmax2^2/ritr_num_Catalyst_Screening);

semilogy(4/3*(1:100:Itr2_SVRG),err_nested_SVRG_hist(1:100:Itr2_SVRG),...
    4/3*(1:100:Itr2_SS),err_nested_SVRG_Screening_hist(1:100:Itr2_SS),...
   ritr_num_Catalyst / Itrmax2 * (1:100:Itr2_Catalyst), err_nested_Catalyst_hist(1:100:Itr2_Catalyst),...
   ritr_num_Catalyst_Screening / Itrmax2 * (1:100:Itr2_CS), err_nested_Catalyst_Screening_hist(1:100:Itr2_CS));
legend('Nested-SVRG','Nested-SVRG-Screening','Nested-Catalyst','Nested-Catalyst-Screening');
xlabel('Number of Component Gradient Evaluations');
ylabel('Relative Error');
ylim([1e-4 10])