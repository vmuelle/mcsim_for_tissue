% lookmcml.m
%	Reads the MCML output file specified (eg., example.mco)
%	and uses the subroutine <getmcml> to read the output file.
%	and return values as the specified globals.
%	If PRINTON ==1 or PLOTON == 1, then getmcml also provides output. 
% USES
%	getmcml.m, which in turn uses 
%		makec2f.m
% by Steven L. Jacques, Jan. 2008
% Oregon Health & Sciences University, Portland, OR, USA

clear global

global A Al Az Azr Fzr Fz Na Nc Nlayers 
global Nr Nz Rd Ra Rr Rra Rsp T Ta Td Tr Tra d dr dz 
global g mua mus n nabove nbelow r z
global rm Azrm Fzrm
global PLOTON PRINTON

PLOTON = 1; PRINTON = 0;

getmcml('data_files/outputs/sample_d_n_470.mco')

