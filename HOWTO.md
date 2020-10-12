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
bundle exec jekyll build
cd _site
mkdir ../../TEMP
cp -r * ../TEMP/
cd ..
git checkout master
cp -r ../TEMP/* .
git add .
git commit -m "Updated the website"
git push
git checkout source
rm -rf ../TEMP
```
