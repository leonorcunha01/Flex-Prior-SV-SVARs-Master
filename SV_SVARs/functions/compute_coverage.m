% function to compute coverage of credible intervals
function [cov90, cov68] = compute_coverage(quant, true_vals)
    cov90 = mean(true_vals >= quant(:,:,1) & true_vals <= quant(:,:,5), 1)';
    cov68 = mean(true_vals >= quant(:,:,2) & true_vals <= quant(:,:,4), 1)';
end
