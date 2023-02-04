library("here")
library("withr")

with_dir(here(), {
    rmarkdown::render("README.Rmd", "html_document")
    system("mv README.html ~/Dropbox/Code/spatialDLPFC_website/index.html")
})

with_dir(
    "~/Dropbox/Code/spatialDLPFC_website",
    system("git ci -am -'Updated website with code/update_website.R'; git push origin gh-pages")
)
