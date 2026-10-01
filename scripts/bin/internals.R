
## Plotting -------------
theme_custom <- function(...) {
  ggpubr::theme_pubr(legend = "right") +
    ggpubr::rotate_x_text(angle = 60) +
    theme(
      plot.title    = element_text(face = "bold"),
      plot.subtitle = element_text(face = "italic"),
      plot.caption  = element_text(face = "italic"),
      strip.text.x    = element_text(face = "bold"),
      strip.text.y    = element_text(face = "bold",angle=0),
      panel.grid    = element_blank(),
      strip.background = element_rect(fill="white",color="grey30")
    )
}



scale_fill_rb <- function(...,limits=NULL) {
  scale_fill_gradient2(
    low      = "#101563",
    mid      = "white",
    high     = "#6b0e20",
    midpoint = 0,
    limits = limits,
    oob  = scales::squish
  )
}


scale_color_rb <- function(...) {
  scale_color_gradient2(
    low      = "#101563",
    mid      = "white",
    high     = "#6b0e20",
    midpoint = 0
  )
}

update_geom_defaults("col",  list(color = "black", linewidth = 0.005))
update_geom_defaults("tile", list(color = "black", linewidth = 0.005))

make_donut <- function(df, fill_var, fill_vals, title) {
  df |>
    dplyr::count(.data[[fill_var]]) |>
    ggplot(aes(
      x    = 2,
      y    = n,
      fill = .data[[fill_var]]
    )) +
    geom_col(
      width     = 1,
      color     = "black",
      linewidth = 0.4
    ) +
    coord_polar(theta = "y") +
    xlim(0.5, 2.5) +
    scale_fill_manual(values = fill_vals) +
    labs(
      #title = title,
      fill  = title
    ) +
    theme_void() +
    theme(
      plot.title      = element_text(hjust = 0.5, size = 11, face = "bold"),
      text = element_text(size=14),
      legend.position = "right"
    )
}

# Might be nice to use ---------

palettes <- list(
  jupiter          = c("#ffd89b", "#19547b"),
  light_bliss      = c("#360033", "#0b8793", "#FFFDE4"),
  mako             = c("#0B0405FF","#357BA2FF","#DEF5E5FF"),
  magma            = c("#000004FF","#B63679FF","#FCFDBFFF"),
  muted_blue       = c("#bbd2c5", "#536976", "#292e49"),
  blue_yellow      = c("#003973", "#e5e5be"),
  lagoon           = c("#43c6ac", "#191654"),
  syracuse         = c("#ffa17f", "#00223e"),
  relay            = c("#3a1c71", "#d76d77", "#ffaf7b"),
  purple_paradise  = c("#1d2b64", "#f8cdda"),
  argon            = c("#03001e", "#7303c0", "#ec38bc", "#fdeff9"),
  celestial        = c("#c33764", "#1d2671"),
  burnt_rose       = c("#eddcd5", "#4a0422"),
  blue_rest        = c("#19547b","grey90","#191654")
)


discrete_palettes <-  list(
  boos_palace = c(
    "#C2697FFF",
    "#FDDFA4FF", 
    "#B4B9E0FF", 
    "#513965FF", 
    "#201B43FF", 
    "#2C3778FF"),
  
  xmas = c(
    "#33645FFF",
    "#9D1B1FFF",
    "#CEA0B5FF", 
    "#F6F5F1FF",
    "#7B7987FF",
    "#1D2731FF"),
  
  dakota = c(
    "#A8CDECFF", 
    "#F6955EFF", 
    "#682C37FF",
    "#9B6981FF",
    "#7887A4FF",
    "#A89F8EFF"),
  
  bayala = c(
    "#F6B6C2FF",
    "#A6A0F2FF",
    "#FBE7B8FF",
    "#6DEDDDFF",
    "#3E6992FF",
    "#FCD47CFF"),
  
  
  kirby = c(
    "#B8396BFF",
    "#FFD1D7FF", 
    "#FFF5CCFF", 
    "#76BAE0FF", 
    "#B28F81FF", 
    "#54483EFF"),
  
  saffron_meadow = c(
    "#E3C1CBFF",
    "#AD5A6BFF", 
    "#C993A2FF", 
    "#365C83FF", 
    "#384351FF", 
    "#4D8F8BFF", 
    "#CDD6ADFF"),
  
  steel_mill = c(
    "#FBE697FF", 
    "#F3AE6DFF", 
    "#516888FF", 
    "#C9DACAFF", 
    "#14232AFF", 
    "#557780FF", 
    "#1F304AFF", 
    "#802729FF"),
  
  
  cyber_punk = c(
    "#B8B8C8FF",
    "#8080A8FF", 
    "#404058FF", 
    "#E0E0E8FF", 
    "#585878FF", 
    "#580800FF", 
    "#F8F8F8FF", 
    "#A83828FF"),
  
  
  periwinkle_dream = c(
    "#D0B8C8FF",
    "#E0C8E0FF", 
    "#284060FF",
    "#F0E0F0FF", 
    "#A090A0FF",
    "#586888FF", 
    "#182848FF", 
    "#686868FF",
    "#90A0D0FF",
    "#F8F8F8FF",
    "#98A0C0FF", 
    "#B0C8F8FF", 
    "#C8E0F8FF"),
  
  akala_forest = c(
    "#D06040FF",
    "#703830FF",
    "#F8D0B0FF",
    "#D8B898FF",
    "#88C0F8FF",
    "#E88058FF",
    "#F8E0D0FF", 
    "#90A8C0FF", 
    "#686078FF",
    "#000000FF",
    "#E08088FF", 
    "#B05040FF", 
    "#A06850FF", 
    "#F8A080FF"),
  
  
  myhero_academia = c(
    "#484850FF", 
    "#000000FF", 
    "#A0B0C0FF", 
    "#604830FF", 
    "#E05060FF",
    "#D8B850FF",
    "#7888A0FF",
    "#F8F8F8FF",
    "#C83850FF", 
    "#888890FF", 
    "#883048FF", 
    "#F8E050FF", 
    "#D8D8D8FF", 
    "#F08070FF", 
    "#485068FF"),
  
  
  pleisosaur = c(
    "#5888C0FF",
    "#3870A8FF", 
    "#80B0E0FF",
    "#000000FF", 
    "#D8E0F8FF",
    "#305078FF",
    "#F8F8F8FF", 
    "#C0C8D8FF", 
    "#F8D880FF", 
    "#A08058FF",
    "#F8F8A0FF", 
    "#F8C060FF"),
  
  
  emperor_penguin =c(
    "#6890D0FF",
    "#E8C8A8FF",
    "#F8E0B8FF",
    "#F8F8F8FF", 
    "#385070FF",
    "#88A8F8FF", 
    "#000000FF",
    "#9898B0FF", 
    "#587890FF", 
    "#D8D8F8FF", 
    "#B8B8D0FF", 
    "#D0B090FF",
    "#B89878FF", 
    "#4068A0FF",
    "#B8C8F8FF"))
