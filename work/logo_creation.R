library(ggplot2)
library(hexSticker)

p <- ggplot(data = data.frame(x = 1:10, y = 1:10), aes(x, y)) + 
  geom_point() + 
  theme_void()

sticker(
  subplot = p, 
  package = "clustGLMM", 
  p_size = 20, 
  s_x = 1, 
  s_y = 0.75, 
  s_width = 1.3, 
  s_height = 1,
  filename = "man/figures/logo_sticker.png"
)
