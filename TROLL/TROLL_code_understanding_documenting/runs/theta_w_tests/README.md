# Hydraulic Test Set – Theta-w Sensitivity

**Creation date:** April 14, 2026

## Purpose of these tests

I was encountering numerical issues where `theta_w` was becoming zero or negative.  
When vegetation was enabled, this led to `NaN` values in photosynthesis and other ecophysiological processes.

These tests were designed to evaluate the model sensitivity to low `theta_w` values and to assess more robust lower-bound constraints for numerical stability.