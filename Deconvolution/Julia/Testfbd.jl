using Conv
using FocusedBlindDecon
using Gadfly

ntg=30 # number of time samples in `g`
nr=20 # number of receivers
nt=40*ntg # time samples in records `d`
nts=nt # samples in `s`

gobs=zeros(ntg,nr) # allocate
#FBD.toy_direct_green!(gobs, c=4.0, bfrac=0.20, afrac=1.0); # add arrival 1
#FBD.toy_reflec_green!(gobs, c=1.5, bfrac=0.35, afrac=-0.6); # add arrival 2
plotg=(x,args...)->spy(x, Guide.xlabel("channel"), Guide.ylabel("time"),args...) # define a plot recipe
p1=plotg(gobs, Guide.title("True g"))

sobs=randn(nts)
plot(y=sobs,x=1:nts,Geom.line, Guide.title("arbitrary source"), Guide.xlabel("time"))

S=Conv.S(sobs, gsize=[ntg,nr], dsize=[nt,nr], slags=[nts-1,0]);
dobs=reshape(S*vec(gobs), nt, nr);

pa=P_fbd(ntg, nt, nr, nts, dobs=dobs, gobs=gobs, sobs=sobs)

FBD.lsbd!(pa)

p2=plotg(pa[:g], Guide.title("LSBD g"))

FBD.fbd!(pa)

p3=plotg(pa[:g], Guide.title("FBD g"))


