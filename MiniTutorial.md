### Если проекта в хранилище ещё совсем нет: ###

1. Создать проект  C# "ClassLibrary".

2. Добавить в проект все исходники с интерфейсами.

3. Проверить что все интерфейсы, которыми будут пользоваться участники других подпроектов, описаны как
"public interface <имя интерфейса>", а не "interface <имя интерфейса>".

4. Проверить, что проект компилируются.

5. Установить какой-нибудь клиент SVN.

Если TortoiseSVN:

5.1 Извлечь рабочую копию: http://tortoisesvn.net/docs/nightly/TortoiseSVN_ru/tsvn-dug-checkout.html.

5.2 Скопировать туда проект (у нас желательно в папку InterfaceProjects).

5.3 Сохранить изменения: http://tortoisesvn.net/docs/nightly/TortoiseSVN_ru/tsvn-dug-commit.html#tsvn-dug-commit-changelists.

### Если проект уже в хранилище: ###

1. Выполнить 5.1 (выше) или обновить рабочую копию: http://tortoisesvn.net/docs/nightly/TortoiseSVN_ru/tsvn-dug-update.html.

2. Внести изменения в код.

3. Выполнить 5.3 (выше).