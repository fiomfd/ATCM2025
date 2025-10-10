% metapost_linearmap1.mp
prologues:=3;
% ruler and compas: a^{1/2} for a=2u>1 
outputformat := "png";
outputtemplate := "%j-%2c.%o";
defaultscale:=2.0;
verbatimtex
%&latex
\documentclass{article}
\begin{document}
etex;
beginfig(1);
u=20;
% R^n
fill (0,0)--(5.19u,3u)--(5.19u,7u)--(0,4u)--cycle withcolor (255,255,0);
draw (0,0)--(5.19u,3u)--(5.19u,7u)--(0,4u)--cycle;
draw (0,2u)--(1.73u,3u);
draw (1.732u,u)--(1.73u,3u);
%R^m
fill (10u,0)--(13.46u,2u)--(13.46u,6u)--(10u,4u)--cycle withcolor 0.9;
draw (10u,0)--(13.46u,2u)--(13.46u,6u)--(10u,4u)--cycle;
% mapping
drawarrow (0,0)--(10u,0) withcolor green;
drawarrow (0,2u)--(10u,2u) withcolor green;
drawarrow (0,4u)--(10u,4u) withcolor green;
drawarrow (1.73u,u)--(10u,0) withcolor green;
drawarrow (1.73u,3u)--(10u,2u) withcolor green;
drawarrow (5.19u,3u)--(10u,0) withcolor green;
drawarrow (5.19u,7u)--(10u,4u) withcolor green;
%
%
pickup pencircle scaled 4pt;
% w-axis
drawarrow (0,0)--(5.19u,3u) withcolor blue;
% v-axis
drawarrow (0,0)--(0,4u) withcolor red;
% z-axis
drawarrow (10u,0)--(13.46u,2u);
% y-axis
drawarrow (10u,0)--(10u,4u) withcolor red+blue;
%
% label
label(btex $R^n \ni u \mapsto Au \in R^m$ etex, (1.5u,6.5u));
label(btex $R^n$ etex, (4.7u,6u));
label(btex $N(A)$ etex, (6u,3.4u));
label(btex $C(A^T)$ etex, (-0.2u,4.6u));
label(btex $R^m$ etex, (13u,5u));
label(btex $N(A^T)$ etex, (14.5u,2u));
label(btex $C(A)$ etex, (10.2u,4.7u));
dotlabel.top(btex $u$ etex, (1.73u,3u));
dotlabel.lft(btex $v$ etex, (0,2u));
dotlabel.bot(btex $w$ etex, (1.73u,u));
dotlabel.bot(btex $0$ etex, (0,0));
dotlabel.rt(btex $Au$ etex, (10u,2u));
dotlabel.bot(btex $0$ etex, (10u,0));
endfig;
end
