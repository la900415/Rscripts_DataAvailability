library(dplyr)
library(ggplot2)
library(tidyr)
################################################################################
################  Mirror bar plot in y-axis   ##################################
################################################################################


## Calling geom_col twice
mpg %>% 
  group_by(manufacturer) %>% 
  summarise_if(is.numeric, mean) %>% 
  select(manufacturer, cty, hwy) %>% 
  ggplot(  ) +
  geom_col(aes(y=manufacturer, x=hwy, fill="hwy")) +
  geom_col(aes(y=manufacturer, x=-cty, fill="cty")) +
  theme_classic( ) +
  theme(axis.text.x = element_text() )

## Converting data to long format
mpg %>% 
  group_by(manufacturer) %>% 
  summarise_if(is.numeric, mean) %>% 
  mutate(cty = -cty) %>% 
  select(manufacturer, cty, hwy) %>% 
  pivot_longer(., -manufacturer) %>% 
  ggplot() +
  geom_col(aes(y=manufacturer, x=value, fill=name)) +
  theme(axis.text.y = element_text(angle = 0, vjust = 0.4))


