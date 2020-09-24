library(ggpubr)

# Regularization
Phi_WT_fixed = function(X, beta, T_mat, W, lambda1=0.2, lambda2=0.01, lambda3=0.1){
  K = ncol(W)
  p = ncol(X)
  N = nrow(X)

  lasso = sum(abs(beta))
  glasso = sum(sqrt((t(W!=0)%*%beta^2)*(apply(t(W!=0), 1, sum))))
  pca = norm(X%*%diag(beta) - T_mat%*%t(W), "F")

  return(lasso*lambda1 + lambda2*glasso + lambda3/2*pca)
}

set.seed(0)
# Generate X
N = 100
X = MASS::mvrnorm(N, mu = c(0, 0, 0),
                  Sigma = matrix(c(1, 0, 0.5, # ... not really important in this example
                                   0, 1, 0,
                                   0.5, 0, 1), nrow = 3, byrow = T))
X <- scale(X)
beta_optim = c(0.5, 0.25, 0.1)

M <- X%*%diag(beta_optim)

svd <- svd(M)
svd$d
T_mat <- svd$u[,1:2]
W <- svd$v[,1:2]%*%diag(svd$d[1:2])
W <- t(apply(W, 1, function(x){x*(abs(x)==max(abs(x)))}))
t(W)%*%W #diagonal
t(T_mat)%*%T_mat # identity

# Create the grid
betas = seq(-1, 1.5, length.out = 100)
grid = expand.grid(beta1 = betas,
                   beta2 = betas, KEEP.OUT.ATTRS = F)
for (i in 1:nrow(grid)) {
  beta = c(grid$beta1[i], grid$beta2[i], 0.1)
  grid$z[i] = Phi_WT_fixed(X, beta, T_mat, W, 2, 1, 1)
}

# plot
gg1 = ggplot(grid) +
  geom_contour(aes(x = beta1, y = beta2, z = z), bins = 12, color='gray')+
  geom_point(aes(x=0.5, y=0.25))+
  geom_point(aes(x=0, y=0))+
  theme_void()

gg1

library(showtext)
## Loading Google fonts (http://www.google.com/fonts)
font_add_google("Share Tech Mono", "Share Tech Mono")
font_add_google("IBM Plex Mono", "ibmplex")
## Automatically use showtext to render text for future devices
showtext_auto()

p <- gg1 + theme_void() + theme_transparent()
sticker(p,
        package="glasp", p_color = "#303030", p_size = 45, p_x = 1, p_y = 1, p_family = "Share Tech Mono",
        s_width = 2, s_height = 2, s_x = 1, s_y = 1,
        h_fill = "#FFFFFF", h_color = "#303030",
        spotlight = F, l_x = 1, l_y = 0.75, l_width = 4, l_height = 1, l_alpha = 0.4,
        url = "github.com/jlaria/glasp", u_size = 5, u_family = "Share Tech Mono")

