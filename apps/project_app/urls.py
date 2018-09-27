from django.conf.urls import url
from . import views
urlpatterns = [
    url(r'^$', views.index),
    url(r'^observe$', views.observe),
    url(r'^display$', views.display),
    url(r'^error$', views.errorpage)
]