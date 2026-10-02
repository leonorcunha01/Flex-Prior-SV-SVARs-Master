% function to compute coverage of intervals for Sig matrix
function [cov90, cov68] = compute_coverage_Sig(quant, true_vals)
    cov90 = mean(true_vals >= squeeze(quant(:,1,:,:,:)) & true_vals <= squeeze(quant(:,5,:,:,:)), 1);
    cov68 = mean(true_vals >= squeeze(quant(:,2,:,:,:)) & true_vals <= squeeze(quant(:,4,:,:,:)), 1);
end
