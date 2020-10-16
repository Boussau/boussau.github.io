---
layout: home
title: Tutorials
subtitle: RevBayes Tutorials
permalink: /tutorials/
levels:
- General structure of the tutorials
- Tutorials for M2 students in Lyon
- Tutorial on dating with relative constraints
levels_id:
- tutorial_structure
---

This list shows some RevBayes tutorials for learning various aspects of RevBayes and Bayesian phylogenetic analysis.
Each one explicitly walks you through model specification and analysis set-up for different phylogenetic methods.
These tutorials have been written for M2 students in Lyon or to provide an example of how to use relative time constraints for dating in RevBayes.

Please see the {% page_ref format %} guide for details about how to read the tutorials.


{% comment %}
{% include keywords.html input=true %}
{% assign keywords = site.empty_array %}
{% endcomment %}

{% assign levels = site.tutorials | sort:"level","last" | group_by:"level" %}
{% for level in levels %}
{% assign i = forloop.index | minus: 1 %}
{% if page.levels[i] %}
<h3>{{ page.levels[i] }}</h3>{:id="{{ page.levels_id[i] }}"}
{% else %}
<h3>Miscellaneous Tutorials</h3>
{% endif %}

{% assign tutorials = level.items | sort:"order","last" %}

<div class="tutorialbox">
{% for tutorial in tutorials %}
{% if tutorial.index %}

{% comment %}{% assign keywords = tutorial.keywords | concat: keywords %}{% endcomment %}

<div class="tutorial {{ tutorial.keywords | join:' '}}">
{% if tutorial.title-old and tutorial.redirect %}
<a class="title" href="https://github.com/revbayes/revbayes_tutorial/raw/master/tutorial_TeX/{{ tutorial.title-old }}/{{ tutorial.title-old }}.pdf">{{ tutorial.title | markdownify }}</a>
{% else %}
<a class="title" href="{{ site.baseurl }}{{ tutorial.url }}">{{ tutorial.title | markdownify }}</a>
{% endif %}
<p class="subtitle" >{{ tutorial.subtitle | markdownify }}</p>
</div>
{% endif %}

{% endfor %}
</div>

{% endfor %}

{% comment %}{% include keywords.html script=true %}{% endcomment %}
