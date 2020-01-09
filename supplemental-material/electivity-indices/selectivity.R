library(selectapref)
library(here)
#download.packages(pkgs="selectapref",destdir = "~/Desktop",type = "source")

ivlev <- function(available, consumed, jacob = FALSE) {
   r <- consumed/sum(consumed)
   p <- available/sum(available)
   if(jacob == TRUE) {
      return ((r-p)/(r+p-(2*r*p)))
   }
   else {
      return ((r-p)/(r+p))
   }
}

eaten <- read.csv(here::here("supplemental-material", "electivity-indices","elect.csv"))
offer <- read.csv(here::here("supplemental-material", "electivity-indices","elect2.csv"))

##Ivlev's
"a"<-ivlev(consumed = eaten$Lionfish.1, available = offer$Lionfish.1, jacob = FALSE)
"b"<-ivlev(consumed = eaten$Lionfish.2, available = offer$Lionfish.1, jacob = FALSE)
"c"<-ivlev(consumed = eaten$Lionfish.3, available = offer$Lionfish.1, jacob = FALSE)
"d"<-ivlev(consumed = eaten$Lionfish.4, available = offer$Lionfish.1, jacob = FALSE)
"e"<-ivlev(consumed = eaten$Lionfish.5, available = offer$Lionfish.1, jacob = FALSE)
"f"<-ivlev(consumed = eaten$Lionfish.6, available = offer$Lionfish.1, jacob = FALSE)
"g"<-ivlev(consumed = eaten$Lionfish.7, available = offer$Lionfish.1, jacob = FALSE)
"h"<-ivlev(consumed = eaten$Lionfish.8, available = offer$Lionfish.1, jacob = FALSE)
"i"<-ivlev(consumed = eaten$Lionfish.9, available = offer$Lionfish.1, jacob = FALSE)
"j"<-ivlev(consumed = eaten$Lionfish.12, available = offer$Lionfish.1, jacob = FALSE)
"k"<-ivlev(consumed = eaten$Lionfish.13, available = offer$Lionfish.1, jacob = FALSE)
"l"<-ivlev(consumed = eaten$Lionfish.16, available = offer$Lionfish.1, jacob = FALSE)
"m"<-ivlev(consumed = eaten$Lionfish.17, available = offer$Lionfish.1, jacob = FALSE)
"n"<-ivlev(consumed = eaten$Lionfish.18, available = offer$Lionfish.1, jacob = FALSE)
"o"<-ivlev(consumed = eaten$Lionfish.20, available = offer$Lionfish.1, jacob = FALSE)
"p"<-ivlev(consumed = eaten$Lionfish.21, available = offer$Lionfish.1, jacob = FALSE)
"q"<-ivlev(consumed = eaten$Lionfish.22, available = offer$Lionfish.1, jacob = FALSE)
"r"<-ivlev(consumed = eaten$Lionfish.25, available = offer$Lionfish.1, jacob = FALSE)
"s"<-ivlev(consumed = eaten$Lionfish.26, available = offer$Lionfish.1, jacob = FALSE)
"t"<-ivlev(consumed = eaten$Lionfish.27, available = offer$Lionfish.1, jacob = FALSE)

iv<-data.frame(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t)
iv2<-data.frame(t(iv))
mean(iv2$X3)



##Jacob's Modified
"a3"<-ivlev(consumed = eaten$Lionfish.1, available = offer$Lionfish.1, jacob = TRUE)
"b3"<-ivlev(consumed = eaten$Lionfish.2, available = offer$Lionfish.1, jacob = TRUE)
"c3"<-ivlev(consumed = eaten$Lionfish.3, available = offer$Lionfish.1, jacob = TRUE)
"d3"<-ivlev(consumed = eaten$Lionfish.4, available = offer$Lionfish.1, jacob = TRUE)
"e3"<-ivlev(consumed = eaten$Lionfish.5, available = offer$Lionfish.1, jacob = TRUE)
"f3"<-ivlev(consumed = eaten$Lionfish.6, available = offer$Lionfish.1, jacob = TRUE)
"g3"<-ivlev(consumed = eaten$Lionfish.7, available = offer$Lionfish.1, jacob = TRUE)
"h3"<-ivlev(consumed = eaten$Lionfish.8, available = offer$Lionfish.1, jacob = TRUE)
"i3"<-ivlev(consumed = eaten$Lionfish.9, available = offer$Lionfish.1, jacob = TRUE)
"j3"<-ivlev(consumed = eaten$Lionfish.12, available = offer$Lionfish.1, jacob = TRUE)
"k3"<-ivlev(consumed = eaten$Lionfish.13, available = offer$Lionfish.1, jacob = TRUE)
"l3"<-ivlev(consumed = eaten$Lionfish.16, available = offer$Lionfish.1, jacob = TRUE)
"m3"<-ivlev(consumed = eaten$Lionfish.17, available = offer$Lionfish.1, jacob = TRUE)
"n3"<-ivlev(consumed = eaten$Lionfish.18, available = offer$Lionfish.1, jacob = TRUE)
"o3"<-ivlev(consumed = eaten$Lionfish.20, available = offer$Lionfish.1, jacob = TRUE)
"p3"<-ivlev(consumed = eaten$Lionfish.21, available = offer$Lionfish.1, jacob = TRUE)
"q3"<-ivlev(consumed = eaten$Lionfish.22, available = offer$Lionfish.1, jacob = TRUE)
"r3"<-ivlev(consumed = eaten$Lionfish.25, available = offer$Lionfish.1, jacob = TRUE)
"s3"<-ivlev(consumed = eaten$Lionfish.26, available = offer$Lionfish.1, jacob = TRUE)
"t3"<-ivlev(consumed = eaten$Lionfish.27, available = offer$Lionfish.1, jacob = TRUE)

jacob<-data.frame(a3,b3,c3,d3,e3,f3,g3,h3,i3,j3,k3,l3,m3,n3,o3,p3,q3,r3,s3,t3)
jacob2<-data.frame(t(jacob))
mean(jacob2$X3)



##Manly's Alpha
"a1"<-manlysalpha(initial = offer$Lionfish.1, consumed = eaten$Lionfish.1)
"b1"<-manlysalpha(initial = offer$Lionfish.2, consumed = eaten$Lionfish.2)
"c1"<-manlysalpha(initial = offer$Lionfish.3, consumed = eaten$Lionfish.3)
"d1"<-manlysalpha(initial = offer$Lionfish.4, consumed = eaten$Lionfish.4)
"e1"<-manlysalpha(initial = offer$Lionfish.5, consumed = eaten$Lionfish.5)
"f1"<-manlysalpha(initial = offer$Lionfish.6, consumed = eaten$Lionfish.6)
"g1"<-manlysalpha(initial = offer$Lionfish.7, consumed = eaten$Lionfish.7)
"h1"<-manlysalpha(initial = offer$Lionfish.8, consumed = eaten$Lionfish.8)
"i1"<-manlysalpha(initial = offer$Lionfish.9, consumed = eaten$Lionfish.9)
"j1"<-manlysalpha(initial = offer$Lionfish.12, consumed = eaten$Lionfish.12)
"k1"<-manlysalpha(initial = offer$Lionfish.13, consumed = eaten$Lionfish.13)
"l1"<-manlysalpha(initial = offer$Lionfish.16, consumed = eaten$Lionfish.16)
"m1"<-manlysalpha(initial = offer$Lionfish.17, consumed = eaten$Lionfish.17)
"n1"<-manlysalpha(initial = offer$Lionfish.18, consumed = eaten$Lionfish.18)
"o1"<-manlysalpha(initial = offer$Lionfish.20, consumed = eaten$Lionfish.20)
"p1"<-manlysalpha(initial = offer$Lionfish.21, consumed = eaten$Lionfish.21)
"q1"<-manlysalpha(initial = offer$Lionfish.22, consumed = eaten$Lionfish.22)
"r1"<-manlysalpha(initial = offer$Lionfish.25, consumed = eaten$Lionfish.25)
"s1"<-manlysalpha(initial = offer$Lionfish.26, consumed = eaten$Lionfish.26)
"t1"<-manlysalpha(initial = offer$Lionfish.27, consumed = eaten$Lionfish.27)

manly<-data.frame(a1,b1,c1,d1,e1,f1,g1,h1,i1,j1,k1,l1,m1,n1,o1,p1,q1,r1,s1,t1)
manly2<-data.frame(t(manly))
mean(manly2$X3)


