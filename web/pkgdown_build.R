library(usethis)
library(pkgdown)

# Generate pkgdown website

# Set up package to use pkgdown
# usethis::use_pkgdown(destdir = 'web/pkgdown')

# This does not need to be rerun unless you want to completely
# reconfigure
# Destination can be changed in _pkgdown.yml, however you also
# need to update .Rbuildignore and .gitignore accordingly

# Build website
pkgdown::build_site()

# When making development changes to the site ensure that 
# development mode is devel in _pkgdown.yml

# Pages can be manually updated after the site has been built
# pkgdown::build_home()
# pkgdown::build_articles() # Includes Getting Started
# pkgdown::build_news()
# pkgdown::build_reference()
# pkgdown::build_tutorials()


