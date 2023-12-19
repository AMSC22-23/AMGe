close all
clear all
%test_multigrid.cpp


A = dlmread("gseidel.txt");
B = dlmread('twolvl.txt');
% C = dlmread('tre.txt');
% D = dlmread('quattro.txt');


% Plot della prima curva
semilogy(A(:,1),A(:,2),'r-')
hold on

semilogy(B(:,1),B(:,2), 'b-')
hold on

semilogy(C(:,1),C(:,2), 'g-')
hold on

semilogy(D(:,1),D(:,2), 'm-')

legend('gs','two','three','four')
grid on
