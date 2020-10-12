# testing the webpage locally:
In the terminal, run:

```
bundle exec jekyll serve
```
Then view the webpage at:

```
http://localhost:4000
```

cf https://docs.github.com/en/free-pro-team@latest/github/working-with-github-pages/testing-your-github-pages-site-locally-with-jekyll

# Building the webpages
```
jekyll build
mv -r _site/* ../TEMP/
git checkout master
mv -r ../TEMP/* .
git add *
git commit -m "Updated the website"
git push
```
