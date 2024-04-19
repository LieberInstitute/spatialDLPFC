library("here")
library("withr")

with_dir(here(), {
    rmarkdown::render("README.Rmd", "html_document")
    system(paste0(
        "mv README.html ",
        file.path(
            dirname(here::here()),
            paste0(basename(here::here()), "_website"),
            "index.html"
        )
    ))
})

with_dir(file.path(dirname(here::here()), paste0(basename(here::here()), "_website")), {
    system("git pull")
    system("git commit -am -'Updated website with code/update_website.R'")
    system("git push origin gh-pages")
})
